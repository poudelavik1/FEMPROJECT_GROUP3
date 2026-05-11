function [ELEMENT, E_tan_elem] = Hex8_ComputeStrainStress(XYZCoord, ELEMCon, ELEMENT, u, E0, nu, sigma_max, varargin)
%HEX8_COMPUTESTRAINSTRESS  Gauss-point strains, stresses, and tangent E
%
%  Called at the end of every Newton-Raphson iteration.
%  Loops over all elements and Gauss points, computes:
%    ε  (6-component strain vector)
%    σ  (6-component stress vector, via tangent D)
%    E_tan (scalar tangent modulus from NonlinearMaterial)
%
%  The tangent E_tan feeds back into the next stiffness assembly.
%
%  Inputs
%    XYZCoord  : [nNode × 3] nodal coordinates
%    ELEMCon   : [NE × 8]   connectivity (1-based local indices)
%    ELEMENT   : structure array (fields added: .strain .stress .E_tan)
%    u         : [3·nNode × 1] current global displacement vector
%    E0        : initial Young's modulus [Pa]
%    nu        : Poisson's ratio
%    sigma_max : saturation stress for NonlinearMaterial [Pa]
%
%  Outputs
%    ELEMENT      : updated structure with per-element Gauss-point data
%    E_tan_elem   : [NE × 1] element-averaged tangent modulus (for assembly)


E_floor = 1e-4 * E0;

NE = size(ELEMCon, 1);
E_tan_elem = zeros(NE, 1);

%% Gauss points and sign tables (same as Hex8_ElementStiffness)
g  = 1/sqrt(3);
gp = [-g -g -g;  g -g -g;  g  g -g; -g  g -g;
      -g -g  g;  g -g  g;  g  g  g; -g  g  g];

xi_i   = [-1  1  1 -1 -1  1  1 -1];
eta_i  = [-1 -1  1  1 -1 -1  1  1];
zeta_i = [-1 -1 -1 -1  1  1  1  1];

for eNo = 1:NE
    conn   = ELEMCon(eNo, :);
    coords = XYZCoord(conn, :);          % [8 × 3]

    %% Build element DOF index vector
    dofIdx = zeros(1, 24);
    for k = 1:8
        n = conn(k);
        dofIdx(3*k-2 : 3*k) = [3*n-2, 3*n-1, 3*n];
    end
    u_e = u(dofIdx);                     % [24 × 1] element displacements

    %% Gauss-point loop
    eps_gp    = zeros(8, 6);
    sig_gp    = zeros(8, 6);
    Etan_gp   = zeros(8, 1);
    sigSAT_gp = zeros(8, 1);   % Saturation equivalent uniaxial stress [Pa]

    for ig = 1:8
        xi   = gp(ig,1);
        eta  = gp(ig,2);
        zeta = gp(ig,3);

        %% Shape function derivatives in natural coordinates
        dN = zeros(3, 8);
        for a = 1:8
            dN(1,a) = 1/8 * xi_i(a)   * (1 + eta_i(a)*eta)   * (1 + zeta_i(a)*zeta);
            dN(2,a) = 1/8 * eta_i(a)  * (1 + xi_i(a)*xi)     * (1 + zeta_i(a)*zeta);
            dN(3,a) = 1/8 * zeta_i(a) * (1 + xi_i(a)*xi)     * (1 + eta_i(a)*eta);
        end

        J     = dN * coords;              % Jacobian [3×3]
        detJ  = det(J);
        if detJ <= 0
            warning('Element %d, GP %d: non-positive Jacobian (%.4g).', eNo, ig, detJ);
        end
        dNxyz = J \ dN;                   % physical derivatives [3×8]

        %% B matrix [6 × 24]
        B = zeros(6, 24);
        for a = 1:8
            c  = 3*(a-1) + 1;
            dx = dNxyz(1,a); dy = dNxyz(2,a); dz = dNxyz(3,a);
            B(:, c:c+2) = [dx 0  0;
                           0  dy 0;
                           0  0  dz;
                           dy dx 0;
                           0  dz dy;
                           dz 0  dx];
        end

        %% Strain vector [εxx εyy εzz γxy γyz γxz]
        eps = B * u_e;                    % [6 × 1]
        eps_gp(ig,:) = eps.';

        %% Effective strain for saturation model
        %  ε_eff = sqrt[ 2/3*(εxx²+εyy²+εzz²-εxx·εyy-εyy·εzz-εxx·εzz)
        %               + 1/2*(γxy²+γyz²+γxz²) ]
        ex = eps(1); ey = eps(2); ez = eps(3);
        gxy= eps(4); gyz= eps(5); gxz= eps(6);

        eps_eff = sqrt( (4/9)*(ex^2 + ey^2 + ez^2 - ex*ey - ey*ez - ex*ez) ...
               + (1/3)*(gxy^2 + gyz^2 + gxz^2) );
        eps_eff = max(eps_eff, 0);

       %% Saturation material: get equivalent stress and tangent E at this Gauss point
[sigma_sat, E_tan] = NonlinearMaterial(eps_eff, E0, sigma_max);

Etan_gp(ig)  = E_tan;
sigSAT_gp(ig) = sigma_sat;

%% Secant modulus for actual stress computation
if eps_eff > 1e-12
    E_sec = sigma_sat / eps_eff;
else
    E_sec = E0;
end
E_sec = max(E_sec, E_floor);
%% Stress tensor using secant D at this Gauss point
D_sec = Hex8_TangentD(E_sec, nu);
sig   = D_sec * eps;                    % [6 × 1]
sig_gp(ig,:) = sig.';
    end

    %% Store results in ELEMENT structure
    ELEMENT(eNo).strain    = eps_gp;     % [8 × 6]
    ELEMENT(eNo).stress    = sig_gp;     % [8 × 6]
    ELEMENT(eNo).E_tan     = Etan_gp;    % [8 × 1]
    ELEMENT(eNo).sigma_sat = sigSAT_gp;  % [8 × 1] saturation equivalent stress [Pa]
    ELEMENT(eNo).sigma_hd  = sigSAT_gp;  % backward-compatible alias for old plotting code
    E_tan_elem(eNo) = mean(abs(Etan_gp));
end

end