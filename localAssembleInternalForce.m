function F_int = localAssembleInternalForce(XYZCoord, ELEMCon, ELEMENT, nNode)

NE   = size(ELEMCon,1);
nDOF = 3*nNode;
F_int = zeros(nDOF,1);

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

for eNo = 1:NE

    conn   = ELEMCon(eNo,:);
    coords = XYZCoord(conn,:);

    dofIdx = zeros(1,24);
    for k = 1:8
        n = conn(k);
        dofIdx(3*k-2:3*k) = [3*n-2, 3*n-1, 3*n];
    end

    fElem = zeros(24,1);

    for ig = 1:8

        xi   = gp(ig,1);
        eta  = gp(ig,2);
        zeta = gp(ig,3);

        dN = zeros(3,8);
        for a = 1:8
            dN(1,a) = 1/8 * xi_i(a)   * (1 + eta_i(a)*eta)   * (1 + zeta_i(a)*zeta);
            dN(2,a) = 1/8 * eta_i(a)  * (1 + xi_i(a)*xi)     * (1 + zeta_i(a)*zeta);
            dN(3,a) = 1/8 * zeta_i(a) * (1 + xi_i(a)*xi)     * (1 + eta_i(a)*eta);
        end

        J     = dN * coords;
        detJ  = det(J);
        dNxyz = J \ dN;

        B = zeros(6,24);
        for a = 1:8
            col = 3*(a-1) + 1;

            dx = dNxyz(1,a);
            dy = dNxyz(2,a);
            dz = dNxyz(3,a);

            B(:,col:col+2) = [dx  0   0;
                0   dy  0;
                0   0   dz;
                dy  dx  0;
                0   dz  dy;
                dz  0   dx];
        end

        sigma = ELEMENT(eNo).stress(ig,:).';

        fElem = fElem + B.' * sigma * detJ;
    end

    F_int(dofIdx) = F_int(dofIdx) + fElem;
end

end