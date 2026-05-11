function ke = Hex8_ElementStiffness(coords, D_input, nu)
ke = zeros(24,24);

g  = 1/sqrt(3);
gp = [-g -g -g;
       g -g -g;
       g  g -g;
      -g  g -g;
      -g -g  g;
       g -g  g;
       g  g  g;
      -g  g  g];

xi_i   = [-1  1  1 -1 -1  1  1 -1];
eta_i  = [-1 -1  1  1 -1 -1  1  1];
zeta_i = [-1 -1 -1 -1  1  1  1  1];

%% Determine D mode: per-GP vector or single matrix
usePerGP = isvector(D_input) && numel(D_input) == 8;

for i = 1:8
    xi   = gp(i,1);
    eta  = gp(i,2);
    zeta = gp(i,3);

    %% Shape function derivatives
    dN = zeros(3,8);
    for a = 1:8
        dN(1,a) = 1/8 * xi_i(a)   * (1 + eta_i(a)*eta)   * (1 + zeta_i(a)*zeta);
        dN(2,a) = 1/8 * eta_i(a)  * (1 + xi_i(a)*xi)     * (1 + zeta_i(a)*zeta);
        dN(3,a) = 1/8 * zeta_i(a) * (1 + xi_i(a)*xi)     * (1 + eta_i(a)*eta);
    end

    J     = dN * coords;
    detJ  = det(J);
    dNxyz = J \ dN;

    %% B matrix [6x24]
    B = zeros(6,24);
    for a = 1:8
        col = 3*(a-1)+1;
        dx = dNxyz(1,a); dy = dNxyz(2,a); dz = dNxyz(3,a);
        B(:,col:col+2) = [dx  0   0;
                           0  dy  0;
                           0   0  dz;
                          dy  dx  0;
                           0  dz  dy;
                          dz   0  dx];
    end

    %% D matrix: per-GP or single
    if usePerGP
        D = Hex8_TangentD(D_input(i), nu);   % GP- D
    else
        D = D_input;                           %  D
    end

    ke = ke + B' * D * B * detJ;
end

end