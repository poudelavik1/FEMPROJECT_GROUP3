function [U_history, R_history, disp_history, ELEMENT, failStep, failReason] = NonlinearFEM_Solver_DisplacementControl_FIXED( ...
    XYZCoord, ELEMCon, NODE, ELEMENT, ...
    E0, nu, sigma_max, ...
    fixedDOF, controlDOF, u_target, ...
    nSteps, maxIter, tol)
%NONLINEARFEM_SOLVER_DISPLACEMENTCONTROL_FIXED
% Outputs:
%   U_history    = displacement vectors, including zero state at column 1
%   R_history    = reaction force at controlled DOF(s), including zero state
%   disp_history = prescribed displacement history, including zero state

if nargin < 12 || isempty(maxIter); maxIter = 50; end       %nargin meacs number of argument input
if nargin < 13 || isempty(tol);     tol     = 1e-6; end

nNode = size(XYZCoord, 1);
NE    = size(ELEMCon,  1);
nDOF  = 3*nNode;

fixedDOF   = unique(fixedDOF(:).');
controlDOF = unique(controlDOF(:).');

% IMPORTANT: controlled DOFs must not also be fixed
if any(ismember(controlDOF, fixedDOF))
    error('A controlDOF is also included in fixedDOF. Check load node/support selection.');
end

knownDOF = unique([fixedDOF, controlDOF]);
freeDOF  = setdiff(1:nDOF, knownDOF);

% Include initial zero state + nSteps nonzero displacement steps
U_history    = zeros(nDOF, nSteps+1);
R_history    = zeros(1,    nSteps+1);
disp_history = zeros(1,    nSteps+1);

u = zeros(nDOF,1);

failStep   = [];
failReason = 'Reached target displacement';

% Initialize element state
for eNo = 1:NE
    ELEMENT(eNo).E_tan    = E0*ones(8,1);
    ELEMENT(eNo).stress   = zeros(8,6);
    ELEMENT(eNo).strain   = zeros(8,6);
    ELEMENT(eNo).sigma_hd = zeros(8,1); % saturation scalar stress alias
end

uSteps = linspace(0, u_target, nSteps+1);

for step = 2:nSteps+1
    targetU = uSteps(step);
    fprintf('\n=== Displacement step %d / %d ===\n', step-1, nSteps);
    fprintf('  Prescribed u = %.6e m = %.3f mm\n', targetU, targetU*1000);

    % Start Newton iteration from previous converged displacement
    u(controlDOF) = targetU;
    u(fixedDOF)   = 0;

    % Update stresses for the initial guess of this step
    [ELEMENT, ~] = Hex8_ComputeStrainStress(XYZCoord, ELEMCon, ELEMENT, u, E0, nu, sigma_max);

    converged = false;

    for iter = 1:maxIter
        % Assemble tangent stiffness using current tangent modulus
        for eNo = 1:NE
            conn   = ELEMCon(eNo,:);
            coords = XYZCoord(conn,:);
            ELEMENT(eNo).stiffness = Hex8_ElementStiffness(coords, ELEMENT(eNo).E_tan, nu);
        end
        K = Hex8_AssembleGlobalStiffness(ELEMENT, ELEMCon, nNode);

        % Assemble current internal force
        F_int = localAssembleInternalForce(XYZCoord, ELEMCon, ELEMENT, nNode);

        % Equilibrium on unknown/free DOFs only: F_ext_free - F_int_free = 0
        % For pure displacement control, F_ext_free = 0.
        Rfree = -F_int(freeDOF);

        resNorm = norm(Rfree);
        refNorm = max(norm(F_int(knownDOF)), 1.0); % reaction scale
        relRes  = resNorm/refNorm;

        fprintf('  Iter %2d  residual = %.3e\n', iter, relRes);

        if relRes < tol
            converged = true;
            fprintf('  Converged in %d iterations\n', iter);
            break
        end

        % Newton correction for free DOFs only
        Kff = K(freeDOF, freeDOF);
        duFree = Kff \ Rfree;

        % Optional safety damping for very large corrections
       maxDu = max(abs(duFree));
        if maxDu > 0.02
            duFree = duFree * (0.02/maxDu);
        end

        u(freeDOF) = u(freeDOF) + duFree;

        % Re-enforce essential BCs exactly
        u(controlDOF) = targetU;
        u(fixedDOF)   = 0;

        % Update stress/tangent for next Newton iteration
        [ELEMENT, ~] = Hex8_ComputeStrainStress(XYZCoord, ELEMCon, ELEMENT, u, E0, nu, sigma_max);
    end

    if ~converged
        warning('Step %d not fully converged. Continuing.',step-1);
    end

    % Final reaction from converged internal force
    F_int = localAssembleInternalForce(XYZCoord, ELEMCon, ELEMENT, nNode);
    reaction = sum(F_int(controlDOF));

    U_history(:,step)    = u;
    R_history(step)      = reaction;
    disp_history(step)   = targetU;

    fprintf('  Reaction = %.3f kN\n', reaction/1000);
end
end

function F_int = localAssembleInternalForce(XYZCoord, ELEMCon, ELEMENT, nNode)
NE   = size(ELEMCon,1);
nDOF = 3*nNode;
F_int = zeros(nDOF,1);

g  = 1/sqrt(3);
gp = [-g -g -g;  g -g -g;  g  g -g; -g  g -g; ...
      -g -g  g;  g -g  g;  g  g  g; -g  g  g];

xi_i   = [-1  1  1 -1 -1  1  1 -1];
eta_i  = [-1 -1  1  1 -1 -1  1  1];
zeta_i = [-1 -1 -1 -1  1  1  1  1];

for eNo = 1:NE
    conn   = ELEMCon(eNo,:);
    coords = XYZCoord(conn,:);

    dofIdx = zeros(1,24);
    for k = 1:8
        n = conn(k);
        dofIdx(3*k-2:3*k) = [3*n-2, 3*n-1, 3*n];
    end

    for ig = 1:8
        xi = gp(ig,1); eta = gp(ig,2); zeta = gp(ig,3);

        dN = zeros(3,8);
        for a = 1:8
            dN(1,a) = 1/8*xi_i(a)  *(1+eta_i(a)*eta)*(1+zeta_i(a)*zeta);
            dN(2,a) = 1/8*eta_i(a) *(1+xi_i(a)*xi)  *(1+zeta_i(a)*zeta);
            dN(3,a) = 1/8*zeta_i(a)*(1+xi_i(a)*xi) *(1+eta_i(a)*eta);
        end

        J     = dN*coords;
        detJ  = det(J);
        dNxyz = J\dN;

        B = zeros(6,24);
        for a = 1:8
            c  = 3*(a-1)+1;
            dx = dNxyz(1,a); dy = dNxyz(2,a); dz = dNxyz(3,a);
            B(:,c:c+2) = [dx 0 0; 0 dy 0; 0 0 dz; dy dx 0; 0 dz dy; dz 0 dx];
        end

        sig = ELEMENT(eNo).stress(ig,:).';
        F_int(dofIdx) = F_int(dofIdx) + B.'*sig*detJ;
    end
end
end
