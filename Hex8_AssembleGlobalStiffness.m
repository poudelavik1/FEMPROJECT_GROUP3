function KG = Hex8_AssembleGlobalStiffness(ELEMENT, ELEMCon, nNode)
%HEX8_ASSEMBLEGLOBALSTIFFNESS Assemble the global stiffness matrix for
% 8-node hexahedral elements with 3 DOF per node: [ux uy uz].
%
% Inputs
%   ELEMENT : structure array with field .stiffness of size [24 x 24]
%   ELEMCon : [NE x 8] connectivity matrix
%   nNode   : total number of nodes
%
% Output
%   KG      : [3*nNode x 3*nNode] global stiffness matrix
%
% Example
%   KG = Hex8_AssembleGlobalStiffness(ELEMENT, ELEMCon, size(XYZCoord,1));

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
