function [XYZCoord, ELEMCon, NODE, ELEMENT] = Hex8_ExtractMeshToExcel(mshFile, excelFile)
%HEX8_EXTRACTMESHTOEXCEL
% Read a Gmsh Hex8 mesh, reorder connectivity as:
%   lower face first, then upper face
% and export to Excel.
%
% lower face = 4 nodes with smaller Z
% upper face = 4 nodes with larger Z
%
% Face ordering inside each face is counterclockwise in the XY plane,
% starting from the "lower-left" point (minimum X+Y).

if nargin < 2 || isempty(excelFile)
    [folder, name, ~] = fileparts(mshFile);
    excelFile = fullfile(folder, [name, '_mesh.xlsx']);
end

[XYZCoord, ELEMCon_raw, elemTags, nodeIDs, elemIDs] = Hex8_ReadGmsh22(mshFile);

nNode = size(XYZCoord,1);
NE    = size(ELEMCon_raw,1);

%% Reorder connectivity: lower face first, upper face second
ELEMCon = zeros(size(ELEMCon_raw));

for eNo = 1:NE
    conn = ELEMCon_raw(eNo,:);
    pts  = XYZCoord(conn,:);   % [8 x 3] = [x y z]

    % --- split by Z: lower 4 nodes and upper 4 nodes
    [~, idxSortZ] = sort(pts(:,3), 'ascend');
    idxLower = idxSortZ(1:4);
    idxUpper = idxSortZ(5:8);

    lowerConn = conn(idxLower);
    upperConn = conn(idxUpper);

    lowerPts = XYZCoord(lowerConn,:);
    upperPts = XYZCoord(upperConn,:);

    % --- order each face counterclockwise in XY
    lowerConn = orderFaceXY(lowerConn, lowerPts);
    upperConn = orderFaceXY(upperConn, upperPts);

    % --- final connectivity
    ELEMCon(eNo,:) = [lowerConn upperConn];
end

%% NODE container
NODE = struct('ID', cell(nNode,1), 'X', [], 'Y', [], 'Z', []);
for i = 1:nNode
    NODE(i).ID = nodeIDs(i);
    NODE(i).X  = XYZCoord(i,1);
    NODE(i).Y  = XYZCoord(i,2);
    NODE(i).Z  = XYZCoord(i,3);
end

%% ELEMENT container
ELEMENT = struct('ID', cell(NE,1), 'con', [], 'X', [], 'Y', [], 'Z', [], 'tags', []);
for eNo = 1:NE
    ELEMENT(eNo).ID   = elemIDs(eNo);
    ELEMENT(eNo).con  = ELEMCon(eNo,:);
    ELEMENT(eNo).X    = XYZCoord(ELEMCon(eNo,:),1).';
    ELEMENT(eNo).Y    = XYZCoord(ELEMCon(eNo,:),2).';
    ELEMENT(eNo).Z    = XYZCoord(ELEMCon(eNo,:),3).';
    ELEMENT(eNo).tags = elemTags(eNo,:);
end

%% Sheet 1: Nodes
nodeTable = table((1:nNode).', nodeIDs, XYZCoord(:,1), XYZCoord(:,2), XYZCoord(:,3), ...
    'VariableNames', {'NodeNo','GmshNodeID','X','Y','Z'});

%% Sheet 2: Elements
varNames = {'ElemNo','GmshElemID','n1','n2','n3','n4','n5','n6','n7','n8'};
elemTable = array2table([(1:NE).', elemIDs, ELEMCon], 'VariableNames', varNames);

if ~isempty(elemTags)
    for k = 1:size(elemTags,2)
        elemTable.(sprintf('tag%d',k)) = elemTags(:,k);
    end
end

%% Sheet 3: Element-node coordinates
rows = NE * 8;
ElemNo = zeros(rows,1);
LocalNode = zeros(rows,1);
GlobalNode = zeros(rows,1);
X = zeros(rows,1);
Y = zeros(rows,1);
Z = zeros(rows,1);

idx = 0;
for eNo = 1:NE
    for a = 1:8
        idx = idx + 1;
        gNode = ELEMCon(eNo,a);
        ElemNo(idx) = eNo;
        LocalNode(idx) = a;
        GlobalNode(idx) = gNode;
        X(idx) = XYZCoord(gNode,1);
        Y(idx) = XYZCoord(gNode,2);
        Z(idx) = XYZCoord(gNode,3);
    end
end

memberCoordTable = table(ElemNo, LocalNode, GlobalNode, X, Y, Z);

%% Write Excel
writetable(nodeTable,        excelFile, 'Sheet', 'Nodes');
writetable(elemTable,        excelFile, 'Sheet', 'Elements');
writetable(memberCoordTable, excelFile, 'Sheet', 'ElementNodeCoords');

fprintf('Mesh data exported to: %s\n', excelFile);
fprintf('Total nodes      : %d\n', nNode);
fprintf('Total Hex8 elems : %d\n', NE);

end


%% --------------------------------------------------------
function connOrdered = orderFaceXY(connFace, ptsFace)
% Order 4 nodes counterclockwise in XY plane
% starting from the point with minimum (X+Y)

xy = ptsFace(:,1:2);
c  = mean(xy,1);

ang = atan2(xy(:,2)-c(2), xy(:,1)-c(1));
[~, order] = sort(ang, 'ascend');   % CCW

connOrdered = connFace(order);
xyOrdered   = xy(order,:);

% rotate so first node is the "lower-left" one
[~, startIdx] = min(xyOrdered(:,1) + xyOrdered(:,2));
connOrdered = circshift(connOrdered, -(startIdx-1));
end