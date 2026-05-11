function KG = Hex8_AssembleGlobalStiffness(ELEMENT, ELEMCon, nNode)
NE   = size(ELEMCon,1);
nDof = nNode * 3;
KG   = zeros(nDof, nDof);

for eNo = 1:NE
    if ~isfield(ELEMENT(eNo), 'stiffness')
        error('ELEMENT(%d) does not contain the field .stiffness', eNo);
    end

    kElem = ELEMENT(eNo).stiffness;

    if ~isequal(size(kElem), [24 24])
        error('ELEMENT(%d).stiffness must be 24x24 for a Hex8 element.', eNo);
    end

    for j = 1:8
        for i = 1:8
            n = ELEMCon(eNo,i);
            m = ELEMCon(eNo,j);

            row = 3*n-2 : 3*n;
            col = 3*m-2 : 3*m;
            kRow = 3*i-2 : 3*i;
            kCol = 3*j-2 : 3*j;

            KG(row,col) = KG(row,col) + kElem(kRow,kCol);
        end
    end
end
end
