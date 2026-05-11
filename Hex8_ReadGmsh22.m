function [XYZCoord, ELEMCon, elemTags, nodeIDs, elemIDs] = Hex8_ReadGmsh22(mshFile)
%HEX8_READGMSH22  Read Hex8 mesh from Gmsh ASCII .msh file.
%                 Supports BOTH format 2.2 and format 4.1 automatically.
%
% Outputs
%   XYZCoord : [nNode x 3] nodal coordinates
%   ELEMCon  : [nElem x 8] connectivity (local 1-based row indices)
%   elemTags : [nElem x nTag] element tags
%   nodeIDs  : [nNode x 1]  original Gmsh node IDs
%   elemIDs  : [nElem x 1]  original Gmsh element IDs

fid = fopen(mshFile, 'r');
if fid == -1
    error('Could not open file: %s', mshFile);
end
cleanup = onCleanup(@() fclose(fid));

%% Detect format version
version = 2.2;   % default
while ~feof(fid)
    line = strtrim(fgetl(fid));
    if strcmp(line, '$MeshFormat')
        fmt_line = strtrim(fgetl(fid));
        version  = sscanf(fmt_line, '%f', 1);
        break;
    end
end
frewind(fid);

fprintf('  Gmsh format detected: %.1f\n', version);

if version >= 4.0
    [XYZCoord, ELEMCon, elemTags, nodeIDs, elemIDs] = readV4(fid);
else
    [XYZCoord, ELEMCon, elemTags, nodeIDs, elemIDs] = readV2(fid);
end

end

%% =========================================================================
%  FORMAT 2.2 READER
%% =========================================================================
function [XYZCoord, ELEMCon, elemTags, nodeIDs, elemIDs] = readV2(fid)

XYZCoord = []; ELEMCon = []; elemTags = []; nodeIDs = []; elemIDs = [];

while ~feof(fid)
    line = strtrim(fgetl(fid));
    if ~ischar(line); break; end

    switch line
        case '$Nodes'
            nNode    = str2double(strtrim(fgetl(fid)));
            nodeIDs  = zeros(nNode, 1);
            XYZCoord = zeros(nNode, 3);
            for i = 1:nNode
                vals        = sscanf(fgetl(fid), '%f')';
                nodeIDs(i)  = vals(1);
                XYZCoord(i,:) = vals(2:4);
            end

        case '$Elements'
            nElemAll   = str2double(strtrim(fgetl(fid)));
            tmpConn    = [];
            tmpTags    = {};
            tmpElemIDs = [];
            for i = 1:nElemAll
                vals     = sscanf(fgetl(fid), '%f')';
                elemType = vals(2);
                nTags    = vals(3);
                if elemType == 5   % Hex8
                    tmpElemIDs(end+1,1) = vals(1);         
                    tmpTags{end+1,1}    = vals(4:3+nTags); 
                    tmpConn(end+1,:)    = vals(4+nTags:end); 
                end
            end
            [ELEMCon, elemTags, elemIDs] = buildConn(tmpConn, tmpTags, tmpElemIDs, nodeIDs);
    end
end
end

%% =========================================================================
%  FORMAT 4.1 READER
%% =========================================================================
function [XYZCoord, ELEMCon, elemTags, nodeIDs, elemIDs] = readV4(fid)

XYZCoord = []; ELEMCon = []; elemTags = []; nodeIDs = []; elemIDs = [];

while ~feof(fid)
    line = strtrim(fgetl(fid));
    if ~ischar(line); break; end

    switch line

        %% --- Nodes ---
        case '$Nodes'
            hdr      = sscanf(fgetl(fid), '%d')';
            % hdr = [numEntityBlocks  numNodes  minNodeTag  maxNodeTag]
            nBlocks  = hdr(1);
            nNode    = hdr(2);
            nodeIDs  = zeros(nNode, 1);
            XYZCoord = zeros(nNode, 3);
            ptr = 0;

            for b = 1:nBlocks
                blk_hdr  = sscanf(fgetl(fid), '%d')';
                % blk_hdr = [entityDim entityTag parametric numNodesInBlock]
                nInBlock = blk_hdr(4);

                % Read node tags first
                tags = zeros(nInBlock, 1);
                for n = 1:nInBlock
                    tags(n) = sscanf(fgetl(fid), '%d', 1);
                end

                % Then read coordinates
                for n = 1:nInBlock
                    ptr = ptr + 1;
                    xyz = sscanf(fgetl(fid), '%f')';
                    nodeIDs(ptr)    = tags(n);
                    XYZCoord(ptr,:) = xyz(1:3);
                end
            end

        %% --- Elements ---
        case '$Elements'
            hdr      = sscanf(fgetl(fid), '%d')';
            % hdr = [numEntityBlocks numElements minElemTag maxElemTag]
            nBlocks  = hdr(1);
            tmpConn    = [];
            tmpTags    = {};
            tmpElemIDs = [];

            for b = 1:nBlocks
                blk_hdr  = sscanf(fgetl(fid), '%d')';
                % blk_hdr = [entityDim entityTag elementType numElementsInBlock]
                elemType = blk_hdr(3);
                nInBlock = blk_hdr(4);

                for e = 1:nInBlock
                    vals = sscanf(fgetl(fid), '%d')';
                    if elemType == 5   % Hex8 only
                        tmpElemIDs(end+1,1) = vals(1);      
                        tmpTags{end+1,1}    = blk_hdr(2);   
                        tmpConn(end+1,:)    = vals(2:9);    
                    end
                end
            end
            [ELEMCon, elemTags, elemIDs] = buildConn(tmpConn, tmpTags, tmpElemIDs, nodeIDs);
    end
end
end

%% =========================================================================
%  SHARED: map raw node IDs -> row indices in XYZCoord
%% =========================================================================
function [ELEMCon, elemTags, elemIDs] = buildConn(tmpConn, tmpTags, tmpElemIDs, nodeIDs)

if isempty(tmpConn)
    error('No 8-node hexahedral elements (type 5) found in mesh file.');
end

[isFound, loc] = ismember(tmpConn, nodeIDs);
if ~all(isFound, 'all')
    error('Some element node IDs not found in node list.');
end
ELEMCon  = loc;
elemIDs  = tmpElemIDs;

maxTagCount = max(cellfun(@numel, tmpTags));
elemTags    = nan(numel(tmpTags), maxTagCount);
for k = 1:numel(tmpTags)
    elemTags(k, 1:numel(tmpTags{k})) = tmpTags{k};
end

end