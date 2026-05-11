function [ELEMENT, E_tan_elem] = Hex8_ComputeStrainStress(XYZCoord, ELEMCon, ELEMENT, u, E0, nu, sigma_max, eps_0, eps_u)

if nargin < 8 || isempty(eps_0);  eps_0 = 0.002;  end   %#ok<NASGU>
if nargin < 9 || isempty(eps_u);  eps_u = 0.0035; end   %#ok<NASGU>

NE         = size(ELEMCon, 1);
E_tan_elem = zeros(NE, 1);

%% Gauss points (2x2x2 full integration)
g  = 1/sqrt(3);
gp = [-g -g -g;  g -g -g;  g  g -g; -g  g -g;
      -g -g  g;  g -g  g;  g  g  g; -g  g  g];

xi_i   = [-1  1  1 -1 -1  1  1 -1];
eta_i  = [-1 -1  1  1 -1 -1  1  1];
zeta_i = [-1 -1 -1 -1  1  1  1  1];

for eNo = 1:NE
    conn   = ELEMCon(eNo, :);
    coords = XYZCoord(conn, :);          % [8 x 3]

    %% Element DOF indices
    dofIdx = zeros(1, 24);
    for k = 1:8
        n = conn(k);
        dofIdx(3*k-2 : 3*k) = [3*n-2, 3*n-1, 3*n];
    end
    u_e = u(dofIdx);                     % [24 x 1]

    %% ---------------------------------------------------------------
    %  PASS 1 — compute strain at all 8 GPs, store B matrices
    %% ---------------------------------------------------------------
    eps_gp = zeros(8, 6);   % strains at each GP
    B_gp   = zeros(6, 24, 8);  % B matrix at each GP

    for ig = 1:8
        xi   = gp(ig,1);
        eta  = gp(ig,2);
        zeta = gp(ig,3);

        dN = zeros(3, 8);
        for a = 1:8
            dN(1,a) = 1/8 * xi_i(a)   * (1 + eta_i(a)*eta)   * (1 + zeta_i(a)*zeta);
            dN(2,a) = 1/8 * eta_i(a)  * (1 + xi_i(a)*xi)     * (1 + zeta_i(a)*zeta);
            dN(3,a) = 1/8 * zeta_i(a) * (1 + xi_i(a)*xi)     * (1 + eta_i(a)*eta);
        end

        J     = dN * coords;
        detJ  = det(J);
        if detJ <= 0
            warning('Element %d, GP %d: non-positive Jacobian (%.4g).', eNo, ig, detJ);
        end
        dNxyz = J \ dN;

        B = zeros(6, 24);
        for a = 1:8
            c  = 3*(a-1) + 1;
            dx = dNxyz(1,a); dy = dNxyz(2,a); dz = dNxyz(3,a);
            B(:, c:c+2) = [dx  0   0;
                            0  dy   0;
                            0   0  dz;
                           dy  dx   0;
                            0  dz  dy;
                           dz   0  dx];
        end

        eps_gp(ig,:)   = (B * u_e).';
        B_gp(:,:,ig)   = B;
    end
    % Maximum absolute strain component at each GP  [8 x 1]
    max_abs_per_gp = max(abs(eps_gp), [], 2);

    % Element-representative strain = mean over all 8 GPs
    eps_eff = mean(max_abs_per_gp);

    % Single NonlinearMaterial call for the whole element
    [sigma_hd_elem, E_tan_e] = NonlinearMaterial(eps_eff, E0, sigma_max);

    % Secant modulus for 3-D stress
    if eps_eff < 1e-15
        E_sec = E0;
    else
        E_sec = sigma_hd_elem / eps_eff;
    end

    % One D matrix for the whole element
    D_sec = Hex8_TangentD(E_sec, nu);

    %% ---------------------------------------------------------------
    %  PASS 2 — compute stresses at all 8 GPs using the element-uniform D
    %% ---------------------------------------------------------------
    sig_gp   = zeros(8, 6);
    Etan_gp  = E_tan_e * ones(8, 1);   % same E_tan at every GP
    sigHD_gp = sigma_hd_elem * ones(8, 1);

    for ig = 1:8
        sig_gp(ig,:) = (D_sec * eps_gp(ig,:).').';
    end

    %% Store results
    ELEMENT(eNo).strain   = eps_gp;    % [8 x 6]
    ELEMENT(eNo).stress   = sig_gp;    % [8 x 6]  [Pa]
    ELEMENT(eNo).E_tan = E_sec * ones(8,1);
    ELEMENT(eNo).sigma_hd = sigHD_gp;  % [8 x 1]  [Pa]  — uniform per element
    E_tan_elem(eNo)       = E_tan_e;
end

end