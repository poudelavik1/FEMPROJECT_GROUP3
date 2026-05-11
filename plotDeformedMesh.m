function plotDeformedMesh(XYZCoord, ELEMCon, u_final, ELEMENT, varargin)
%PLOTDEFORMEDMESH  Deformed shape with concrete-appropriate stress contours.
%
%  All fields are physically meaningful for concrete.
%  Von Mises is intentionally excluded (it is a steel criterion).
%
% -------------------------------------------------------------------------
%  CONCRETE FIELDS
% -------------------------------------------------------------------------
%  'sigma_hd'       Hognestad uniaxial equivalent stress [MPa]
%                   Same quantity used by the solver for failure.
%                   Limit = f'c (sigma_max).
%
%  'drucker_prager' Drucker-Prager demand/capacity ratio [-]   0 → 1
%                   Proper 3-D pressure-dependent concrete criterion.
%                   Accounts for concrete being ~10x weaker in tension.
%                   D/C = 1.0 means the point is on the failure surface.
%                   Needs 'fc' and 'ft' parameters.
%
%  'principal_max'  Maximum principal stress σ1 [MPa]
%                   Positive = tension → concrete cracks when σ1 > f't.
%                   Red zones = cracking risk.
%
%  'principal_min'  Minimum principal stress σ3 [MPa]
%                   Negative = compression → crushing when σ3 < -f'c.
%                   Blue zones = high compression.
%
%  'pressure'       Hydrostatic pressure p = -I1/3 [MPa]
%                   Positive = net compression (good for concrete).
%                   Negative = net tension (bad for concrete).
%
%  'E_tan'          Tangent modulus [GPa]
%                   Low E_tan = softened/damaged zone.
%
%  'dispMag'        Displacement magnitude [mm]
%  'dispX/Y/Z'      Component displacements [mm]
%
% -------------------------------------------------------------------------
%  USAGE
% -------------------------------------------------------------------------
%  plotDeformedMesh(XYZCoord, ELEMCon, u_final, ELEMENT)
%  plotDeformedMesh(..., 'Field','drucker_prager','fc',30e6,'ft',3e6)
%
%  Name-value options
%    'Field'       string  — field name (default: 'sigma_hd')
%    'fc'          [Pa]    — compressive strength (default: 30e6)
%    'ft'          [Pa]    — tensile strength     (default: 0.1*fc)
%    'ScaleFactor' scalar  — deformation magnification ([] = auto)
%    'nContours'   int     — colour bands (default: 20)
%    'ShowGhost'   logical — undeformed outline (default: true)
%    'Title'       string  — override figure title (default: [])

%% ========================================================================
%  Parse inputs
%% ========================================================================
p = inputParser;
addParameter(p, 'Field',       'sigma_hd');
addParameter(p, 'fc',          30e6);
addParameter(p, 'ft',          []);
addParameter(p, 'ScaleFactor', []);
addParameter(p, 'nContours',   20);
addParameter(p, 'ShowGhost',   true);
addParameter(p, 'Title',       []);
parse(p, varargin{:});

fieldName   = p.Results.Field;
fc          = p.Results.fc;
ft          = p.Results.ft;
scaleFactor = p.Results.ScaleFactor;
nContours   = p.Results.nContours;
showGhost   = p.Results.ShowGhost;
userTitle   = p.Results.Title;

if isempty(ft)
    ft = 0.1 * fc;   % typical: f't ≈ 10% of f'c
end

nNode = size(XYZCoord, 1);
NE    = size(ELEMCon,  1);

%% ========================================================================
%  Displacements
%% ========================================================================
UX = u_final(1:3:end);
UY = u_final(2:3:end);
UZ = u_final(3:3:end);

%% ========================================================================
%  Auto scale factor  — deformation spans ~10% of bounding box
%% ========================================================================
if isempty(scaleFactor)
    maxDisp = max(sqrt(UX.^2 + UY.^2 + UZ.^2));
    maxDim  = max(max(XYZCoord,[],1) - min(XYZCoord,[],1));
    if maxDisp > 1e-20
        scaleFactor = 0.10 * maxDim / maxDisp;
    else
        scaleFactor = 1;
    end
    fprintf('[plotDeformedMesh]  Auto scale factor = %.2f\n', scaleFactor);
end

%% ========================================================================
%  Deformed coordinates
%% ========================================================================
XYZ_def = XYZCoord + scaleFactor * [UX, UY, UZ];

%% ========================================================================
%  Compute nodal field
%% ========================================================================
[nodalField, cLabel, cmapFcn, diverging] = computeNodalField( ...
    fieldName, ELEMENT, UX, UY, UZ, ELEMCon, nNode, NE, fc, ft);

%% ========================================================================
%  Hex8 face definitions  (lower face 1-4 CCW, upper face 5-8 CCW)
%% ========================================================================
hexFaces = [1 2 3 4;   % bottom
            5 6 7 8;   % top
            1 2 6 5;   % front
            2 3 7 6;   % right
            3 4 8 7;   % back
            4 1 5 8];  % left

%% ========================================================================
%  Build global face list
%% ========================================================================
allFaces = zeros(NE*6, 4, 'uint32');
fIdx = 0;
for eNo = 1:NE
    conn = ELEMCon(eNo,:);
    for f = 1:6
        fIdx = fIdx + 1;
        allFaces(fIdx,:) = conn(hexFaces(f,:));
    end
end

%% ========================================================================
%  Figure
%% ========================================================================
fig = figure('Color','w','Name',['Concrete — ' fieldName]);
fig.Position = [80 80 1100 720];

%% Ghost outline (undeformed)
if showGhost
    patch('Vertices',  XYZCoord, ...
          'Faces',     allFaces, ...
          'FaceColor', 'none', ...
          'EdgeColor', [0.78 0.78 0.78], ...
          'LineStyle', '--', ...
          'LineWidth', 0.4, ...
          'FaceAlpha', 0);
    hold on;
end

%% Deformed coloured mesh
hp = patch('Vertices',        XYZ_def, ...
           'Faces',           allFaces, ...
           'FaceVertexCData', nodalField, ...
           'FaceColor',       'interp', ...
           'EdgeColor',       [0.18 0.18 0.18], ...
           'LineWidth',       0.5, ...
           'FaceAlpha',       0.92);

%% Colormap
colormap(cmapFcn(nContours));

%% Colorbar
cb              = colorbar('eastoutside');
cb.Label.String = cLabel;
cb.Label.FontSize = 11;
cb.FontSize     = 10;

%% Color limits
if strcmpi(fieldName, 'drucker_prager')
    clim([0 1]);
    text(0.02, 0.97, 'D/C = 1.0  \rightarrow  Failure surface', ...
        'Units','normalized','FontSize',9,'Color',[0.8 0 0], ...
        'FontWeight','bold','VerticalAlignment','top');
elseif diverging
    absMax = max(abs(nodalField));
    if absMax < 1e-30; absMax = 1; end
    clim([-absMax, absMax]);
else
    cLim = [min(nodalField), max(nodalField)];
    if diff(cLim) < 1e-30; cLim = cLim + [-1 1]; end
    clim(cLim);
end

%% Axes — pbaspect keeps geometry correct without inflating empty space
%  axis equal in 3D expands all axes to the LARGEST span, creating a huge
%  empty cube. pbaspect instead scales the plot box to match data proportions.
allPts = [XYZCoord; XYZ_def];
lo     = min(allPts,[],1);
hi     = max(allPts,[],1);
span   = hi - lo;
pad    = 0.08 * max(span);
lo     = lo - pad;
hi     = hi + pad;
span   = hi - lo;          % recompute after padding

view(3); grid on;
xlim([lo(1), hi(1)]);
ylim([lo(2), hi(2)]);
zlim([lo(3), hi(3)]);
pbaspect([span(1), span(2), span(3)]);   % visual box matches data shape
xlabel('X [m]','FontSize',11);
ylabel('Y [m]','FontSize',11);
zlabel('Z [m]','FontSize',11);
set(gca,'FontSize',10,'Box','on');

if isempty(userTitle)
    titleStr = sprintf('Deformed Shape (\\times%.0f)  |  %s', ...
                       scaleFactor, strrep(fieldName,'_','\_'));
else
    titleStr = userTitle;
end
title(titleStr,'FontSize',13,'FontWeight','bold');

%% Legend
if showGhost
    hold off;
    hGhost = patch('Visible','off','FaceColor','none', ...
                   'EdgeColor',[0.78 0.78 0.78],'LineStyle','--');
    legend([hp, hGhost], ...
           {sprintf('Deformed (\\times%.0f)', scaleFactor), 'Undeformed'}, ...
           'Location','best','FontSize',9);
end

drawnow;
end


%% =========================================================================
%  LOCAL: compute nodal field values
%% =========================================================================
function [nodalField, cLabel, cmapFcn, diverging] = computeNodalField( ...
    fieldName, ELEMENT, UX, UY, UZ, ELEMCon, nNode, NE, fc, ft)

diverging = false;

switch lower(fieldName)

    %----------------------------------------------------------------------
    % Hognestad stress — same as solver failure criterion
    %----------------------------------------------------------------------
    case 'sigma_hd'
        elemVal    = arrayfun(@(e) mean(abs(ELEMENT(e).sigma_hd)), 1:NE).';
        nodalField = elemToNode(elemVal, ELEMCon, nNode) * 1e-6;
        cLabel     = 'Hognestad Stress |σ_{hd}| [MPa]';
        cmapFcn    = @jet;

    %----------------------------------------------------------------------
    % Drucker-Prager D/C ratio — proper 3-D concrete failure criterion
    %
    % Calibrated to uniaxial f'c and f't:
    %   α  = (fc - ft) / ( √3·(fc + ft) )
    %   k  = ft · (α + 1/√3)
    %   DP = α·I1 + √J2
    %   DC = DP / k   (clamped to [0,1]; value = 1 means on failure surface)
    %----------------------------------------------------------------------
    case 'drucker_prager'
        alpha_dp = (fc - ft) / (sqrt(3) * (fc + ft));
        k_dp     = ft * (alpha_dp + 1/sqrt(3));

        elemVal = zeros(NE,1);
        for eNo = 1:NE
            s  = ELEMENT(eNo).stress;           % [8 x 6]
            I1 = s(:,1) + s(:,2) + s(:,3);
            J2 = (1/6)*((s(:,1)-s(:,2)).^2 + ...
                        (s(:,2)-s(:,3)).^2 + ...
                        (s(:,3)-s(:,1)).^2) + ...
                 s(:,4).^2 + s(:,5).^2 + s(:,6).^2;
            dp_stress       = alpha_dp * I1 + sqrt(max(J2, 0));
            dc              = max(dp_stress / k_dp, 0);   % <0 means inside cap
            elemVal(eNo)    = mean(dc);
        end
        nodalField = min(elemToNode(elemVal, ELEMCon, nNode), 1.0);
        cLabel     = sprintf('Drucker-Prager D/C  (f''c = %.0f MPa, f''t = %.0f MPa)', ...
                             fc*1e-6, ft*1e-6);
        cmapFcn    = @jet;

    %----------------------------------------------------------------------
    % Maximum principal stress σ1  (tension indicator)
    % Positive = tension, red = cracking risk
    % Concrete cracks when σ1 > f't
    %----------------------------------------------------------------------
    case 'principal_max'
        elemVal = zeros(NE,1);
        for eNo = 1:NE
            s = ELEMENT(eNo).stress;
            p1 = zeros(8,1);
            for ig = 1:8
                eigv    = eig(stressTensor(s(ig,:)));
                p1(ig)  = max(eigv);
            end
            elemVal(eNo) = mean(p1);
        end
        nodalField = elemToNode(elemVal, ELEMCon, nNode) * 1e-6;
        cLabel     = sprintf('Max Principal Stress σ_1 [MPa]  (red = tension, f''t = %.0f MPa)', ...
                             ft*1e-6);
        cmapFcn    = @coolwarm;
        diverging  = true;

    %----------------------------------------------------------------------
    % Minimum principal stress σ3  (compression indicator)
    % Negative = compression, blue = high compression
    % Concrete crushes when σ3 < -f'c
    %----------------------------------------------------------------------
    case 'principal_min'
        elemVal = zeros(NE,1);
        for eNo = 1:NE
            s = ELEMENT(eNo).stress;
            p3 = zeros(8,1);
            for ig = 1:8
                eigv    = eig(stressTensor(s(ig,:)));
                p3(ig)  = min(eigv);
            end
            elemVal(eNo) = mean(p3);
        end
        nodalField = elemToNode(elemVal, ELEMCon, nNode) * 1e-6;
        cLabel     = sprintf('Min Principal Stress σ_3 [MPa]  (blue = compression, f''c = %.0f MPa)', ...
                             fc*1e-6);
        cmapFcn    = @coolwarm;
        diverging  = true;

    %----------------------------------------------------------------------
    % Hydrostatic pressure  p = -I1/3
    % Positive = net compression (strengthens concrete via confinement)
    % Negative = net tension (weakens concrete)
    %----------------------------------------------------------------------
    case 'pressure'
        elemVal = zeros(NE,1);
        for eNo = 1:NE
            s  = ELEMENT(eNo).stress;
            I1 = s(:,1) + s(:,2) + s(:,3);
            elemVal(eNo) = mean(-I1 / 3);
        end
        nodalField = elemToNode(elemVal, ELEMCon, nNode) * 1e-6;
        cLabel     = 'Hydrostatic Pressure p = -I_1/3 [MPa]  (red = compression +)';
        cmapFcn    = @coolwarm;
        diverging  = true;

    %----------------------------------------------------------------------
    % Tangent modulus — softening/damage indicator
    % Low E_tan (blue end) = material has softened significantly
    % High E_tan (red end) = still near-elastic
    %----------------------------------------------------------------------
    case 'e_tan'
        elemVal    = arrayfun(@(e) mean(ELEMENT(e).E_tan), 1:NE).';
        nodalField = elemToNode(elemVal, ELEMCon, nNode) * 1e-9;
        cLabel     = 'Tangent Modulus E_{tan} [GPa]  (blue = softened / damaged)';
        cmapFcn    = @(n) flipud(jet(n));   % reversed: blue = low E

    %----------------------------------------------------------------------
    % Displacements
    %----------------------------------------------------------------------
    case 'dispmag'
        nodalField = sqrt(UX.^2 + UY.^2 + UZ.^2) * 1e3;
        cLabel     = 'Displacement Magnitude [mm]';
        cmapFcn    = @jet;

    case 'dispx'
        nodalField = UX * 1e3;
        cLabel     = 'X Displacement [mm]';
        cmapFcn    = @coolwarm;
        diverging  = true;

    case 'dispy'
        nodalField = UY * 1e3;
        cLabel     = 'Y Displacement [mm]';
        cmapFcn    = @coolwarm;
        diverging  = true;

    case 'dispz'
        nodalField = UZ * 1e3;
        cLabel     = 'Z Displacement [mm]';
        cmapFcn    = @coolwarm;
        diverging  = true;

    otherwise
        error(['plotDeformedMesh: unknown Field ''%s''.\n' ...
               'Valid: sigma_hd, drucker_prager, principal_max, principal_min,\n' ...
               '       pressure, E_tan, dispMag, dispX, dispY, dispZ'], fieldName);
end
end


%% =========================================================================
%  LOCAL: build 3x3 stress tensor from 6-component stress vector
%  Input order: [sxx syy szz txy tyz txz]
%% =========================================================================
function T = stressTensor(s)
T = [s(1) s(4) s(6);
     s(4) s(2) s(5);
     s(6) s(5) s(3)];
end


%% =========================================================================
%  LOCAL: scatter element values to nodes and average
%% =========================================================================
function nodalVal = elemToNode(elemVal, ELEMCon, nNode)
nodeSum   = zeros(nNode,1);
nodeCount = zeros(nNode,1);
NE = size(ELEMCon,1);
for eNo = 1:NE
    conn = ELEMCon(eNo,:);
    nodeSum(conn)   = nodeSum(conn)   + elemVal(eNo);
    nodeCount(conn) = nodeCount(conn) + 1;
end
nodeCount(nodeCount == 0) = 1;
nodalVal = nodeSum ./ nodeCount;
end


%% =========================================================================
%  LOCAL: blue-white-red diverging colormap  (no toolbox needed)
%  Blue = negative (compression/inward), White = zero, Red = positive (tension)
%% =========================================================================
function cmap = coolwarm(n)
if nargin < 1; n = 256; end
t    = linspace(0, 1, n).';
r    = min(2*t,       1.0);
g    = 1 - abs(2*t - 1);
b    = min(2*(1-t),   1.0);
cmap = [r, g, b];
end