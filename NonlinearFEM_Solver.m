function [U_history, F_history, ELEMENT, failStep, failReason] = NonlinearFEM_Solver( ...
    XYZCoord, ELEMCon, NODE, ELEMENT, ...
    E0, nu, sigma_max, ...
    F_total, fixedDOF, ...
    Fmax_total, nSteps, maxIter, tol, varargin)

%% -----------------------------------------------------------------------
if nargin < 12; maxIter = 100;   end  %#ok<NASGU>
if nargin < 13; tol     = 1e-4; end  %#ok<NASGU>

dispLimit      = 0.058;   % 58 mm
dispMonitorDOF = [];

for k = 1:2:length(varargin)
    switch lower(varargin{k})
        case 'displimit';       dispLimit      = varargin{k+1};
        case 'dispmonitordof';  dispMonitorDOF = varargin{k+1};
        otherwise
            warning('NonlinearFEM_Solver: unknown option ''%s''', varargin{k});
    end
end

if isempty(dispMonitorDOF)
    [~, dispMonitorDOF] = max(abs(F_total));
end

%% Mesh
nNode = size(XYZCoord, 1);
NE    = size(ELEMCon,  1);
nDOF  = 3 * nNode;

allDOF  = (1:nDOF).';
freeDOF = setdiff(allDOF, fixedDOF(:));

U_history = zeros(nDOF, nSteps);
F_history  = zeros(1,   nSteps);
u          = zeros(nDOF, 1);

failStep   = [];
failReason = 'All steps completed or displacement limit reached';

%% Constant incremental load — same dF every step
dF = F_total / nSteps;

%% Initialise element fields
for eNo = 1:NE
    ELEMENT(eNo).E_tan    = E0 * ones(8,1);
    ELEMENT(eNo).stress   = zeros(8,6);
    ELEMENT(eNo).strain   = zeros(8,6);
    ELEMENT(eNo).sigma_hd = zeros(8,1);
end

fprintf('\n[NonlinearFEM_Solver]  PURE INCREMENTAL (Forward-Euler)\n');
fprintf('  Steps      : %d\n',     nSteps);
fprintf('  dF/step    : %.3f kN\n', norm(dF)/1e3);
fprintf('  Disp limit : %.1f mm  (DOF %d)\n\n', dispLimit*1e3, dispMonitorDOF);

%% =========================================================================
%  Main incremental loop
%% =========================================================================
for step = 1:nSteps

    pct = 100 * step / nSteps;
    fprintf('--- Step %d/%d  (%.1f%%) ---\n', step, nSteps, pct);

    %% 1. Build tangent stiffness from current E_tan
    for eNo = 1:NE
        conn   = ELEMCon(eNo, :);
        coords = XYZCoord(conn, :);
        ELEMENT(eNo).stiffness = Hex8_ElementStiffness( ...
            coords, ELEMENT(eNo).E_tan, nu);
    end
    K      = Hex8_AssembleGlobalStiffness(ELEMENT, ELEMCon, nNode);
    K_free = K(freeDOF, freeDOF);

    %% 2. Solve — no rcond abort; let MATLAB's backslash handle conditioning
    du_free = K_free \ dF(freeDOF);

    %% 3. Guard: stop only if solution is numerically broken (NaN / Inf)
    if any(~isfinite(du_free))
        failStep   = step;
        failReason = sprintf( ...
            'Non-finite displacement increment at step %d (%.1f%%) — K truly singular.', ...
            step, pct);
        fprintf('\n*** NON-FINITE du at step %d — structural collapse ***\n%s\n', ...
                step, failReason);
        if step > 1
            U_history = U_history(:, 1:step-1);
            F_history = F_history(:, 1:step-1);
        else
            U_history = zeros(nDOF, 1);
            F_history = 0;
        end
        return
    end

    %% 4. Update displacements
    du          = zeros(nDOF, 1);
    du(freeDOF) = du_free;
    u           = u + du;

    %% 5. Update material state (strain / stress / E_tan / sigma_hd)
    [ELEMENT, E_tan_elem] = Hex8_ComputeStrainStress( ...
        XYZCoord, ELEMCon, ELEMENT, u, E0, nu, sigma_max);

    %% 6. Diagnostics
    maxHD  = max(arrayfun(@(e) max(ELEMENT(e).sigma_hd), 1:NE));
    u_mon  = abs(u(dispMonitorDOF));

    fprintf('  |du|=%.3e m | u_mon=%.3f mm | HD_max=%.3f MPa | E_tan: min=%.4f mean=%.4f GPa\n', ...
            norm(du), u_mon*1e3, maxHD/1e6, min(E_tan_elem)/1e9, mean(E_tan_elem)/1e9);

    %% 7. Store
    U_history(:, step) = u;
    F_history(step)    = step * norm(dF);

    %% 8. Displacement stop criterion
    if u_mon >= dispLimit
        failStep   = step;
        failReason = sprintf( ...
            'Displacement limit %.1f mm reached at step %d (%.1f%% load).  u_mon = %.3f mm.', ...
            dispLimit*1e3, step, pct, u_mon*1e3);
        fprintf('\n*** DISP LIMIT %.1f mm reached at step %d ***\n%s\n', ...
                dispLimit*1e3, step, failReason);
        U_history = U_history(:, 1:step);
        F_history = F_history(:, 1:step);
        return
    end

end

fprintf('\nAll %d steps done.  Final monitor disp = %.3f mm\n', ...
        nSteps, abs(u(dispMonitorDOF))*1e3);
end