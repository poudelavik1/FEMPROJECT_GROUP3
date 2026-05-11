function plotDeformedShape(XYZCoord, ELEMCon, u_final, stressValue, varargin)
%PLOTDEFORMEDSHAPE Plot undeformed + deformed HEX8 mesh with optional stress contour.
%
% Usage:
%   plotDeformedShape(XYZCoord, ELEMCon, u_final)
%   plotDeformedShape(XYZCoord, ELEMCon, u_final, stressValue)
%   plotDeformedShape(..., 'ScaleFactor', 1, 'StressLabel', 'Stress [MPa]')
%
% stressValue can be either:
%   nNode x 1 nodal stress values, or
%   nElem x 1 element stress values.

%% Optional stress input
if nargin < 4
    stressValue = [];
end

%% Parse inputs
p = inputParser;
addParameter(p, 'ScaleFactor', []);
addParameter(p, 'Padding',     0.08);
addParameter(p, 'Title',       []);
addParameter(p, 'StressLabel', 'Stress');
addParameter(p, 'ShowStress',  true);
addParameter(p, 'StressUnit',  'MPa');   % 'Pa' or 'MPa'
parse(p, varargin{:});

scaleFactor = p.Results.ScaleFactor;
padding     = p.Results.Padding;
userTitle   = p.Results.Title;
stressLabel = p.Results.StressLabel;
showStress  = p.Results.ShowStress;
stressUnit  = p.Results.StressUnit;

%% Basic checks
nNode = size(XYZCoord,1);
NE    = size(ELEMCon,1);

if size(XYZCoord,2) ~= 3
    error('XYZCoord must be nNode x 3.');
end
if size(ELEMCon,2) ~= 8
    error('ELEMCon must be nElem x 8 for HEX8 elements.');
end
if numel(u_final) ~= 3*nNode
    error('u_final must have 3*nNode entries.');
end

u_final = u_final(:);

%% Displacements
UX = u_final(1:3:end);
UY = u_final(2:3:end);
UZ = u_final(3:3:end);

%% Auto scale factor
if isempty(scaleFactor)
    maxDisp = max(sqrt(UX.^2 + UY.^2 + UZ.^2));
    maxDim  = max(max(XYZCoord,[],1) - min(XYZCoord,[],1));
    if maxDisp > 1e-20
        scaleFactor = 0.10 * maxDim / maxDisp;
    else
        scaleFactor = 1;
    end
    fprintf('[plotDeformedShape] Auto scale factor = %.2f\n', scaleFactor);
end

%% Deformed coordinates
XYZ_def = XYZCoord + scaleFactor * [UX, UY, UZ];

%% HEX8 faces
hexFaces = [1 2 3 4;
            5 6 7 8;
            1 2 6 5;
            2 3 7 6;
            3 4 8 7;
            4 1 5 8];

allFaces = zeros(NE*6,4);
fIdx = 0;
for eNo = 1:NE
    conn = ELEMCon(eNo,:);
    for f = 1:6
        fIdx = fIdx + 1;
        allFaces(fIdx,:) = conn(hexFaces(f,:));
    end
end

%% Convert stress to nodal values
hasStress = showStress && ~isempty(stressValue);

if hasStress
    stressValue = stressValue(:);

    if numel(stressValue) == NE
        nodalStress = zeros(nNode,1);
        count = zeros(nNode,1);
        for eNo = 1:NE
            conn = ELEMCon(eNo,:);
            nodalStress(conn) = nodalStress(conn) + stressValue(eNo);
            count(conn) = count(conn) + 1;
        end
        nodalStress = nodalStress ./ max(count,1);
    elseif numel(stressValue) == nNode
        nodalStress = stressValue;
    else
        error('stressValue must be either nElem x 1 or nNode x 1.');
    end

    % Convert to MPa only for display if values are in Pa
    switch lower(stressUnit)
        case 'mpa'
            nodalStressPlot = nodalStress / 1e6;
            unitText = 'MPa';
        otherwise
            nodalStressPlot = nodalStress;
            unitText = 'Pa';
    end

    stressMin = min(nodalStressPlot);
    stressMax = max(nodalStressPlot);
else
    nodalStressPlot = [];
    stressMin = [];
    stressMax = [];
    unitText = '';
end

%% Axis limits
allPts = [XYZCoord; XYZ_def];
lo = min(allPts,[],1);
hi = max(allPts,[],1);
span = hi - lo;
pad = padding * max(span);
lo = lo - pad;
hi = hi + pad;
span = hi - lo;

%% Figure
figure('Color','w','Name','Deformed Shape');
set(gcf,'Position',[100 100 1100 720]);
hold on;

%% Undeformed mesh
h_und = patch('Vertices',XYZCoord, ...
              'Faces',allFaces, ...
              'FaceColor',[0.80 0.80 0.80], ...
              'FaceAlpha',0.18, ...
              'EdgeColor',[0.50 0.50 0.50], ...
              'LineStyle','--', ...
              'LineWidth',0.6);

%% Deformed mesh
if hasStress
    h_def = patch('Vertices',XYZ_def, ...
                  'Faces',allFaces, ...
                  'FaceVertexCData',nodalStressPlot, ...
                  'FaceColor','interp', ...
                  'FaceAlpha',0.95, ...
                  'EdgeColor',[0.05 0.05 0.05], ...
                  'LineStyle','-', ...
                  'LineWidth',0.45);
else
    h_def = patch('Vertices',XYZ_def, ...
                  'Faces',allFaces, ...
                  'FaceColor',[0.18 0.44 0.78], ...
                  'FaceAlpha',0.75, ...
                  'EdgeColor',[0.05 0.20 0.45], ...
                  'LineStyle','-', ...
                  'LineWidth',0.8);
end

%% Axes
view(3); grid on; axis equal;
xlim([lo(1),hi(1)]);
ylim([lo(2),hi(2)]);
zlim([lo(3),hi(3)]);
pbaspect([span(1),span(2),span(3)]);
xlabel('X [m]','FontSize',11);
ylabel('Y [m]','FontSize',11);
zlabel('Z [m]','FontSize',11);
set(gca,'FontSize',10,'Box','on');

%% Colorbar and visible max/min annotation
if hasStress
    colormap(jet(256));
    if abs(stressMax - stressMin) > eps
        caxis([stressMin stressMax]);
    end

    cb = colorbar;
    cb.Label.String = sprintf('%s [%s]', stressLabel, unitText);
    cb.Label.FontSize = 11;

    % Display max and min stress on the figure
    maxStress = stressMax;
    minStress = stressMin;

   annotationText = sprintf('Max Stress = %.4f MPa\nMin Stress = %.4f MPa', ...
                         stressMax/1e6, stressMin/1e6);

    annotation('textbox', [0.72 0.78 0.20 0.10], ...
        'String', annotationText, ...
        'FitBoxToText', 'on', ...
        'BackgroundColor', 'w', ...
        'EdgeColor', 'k', ...
        'FontSize', 10, ...
        'FontWeight', 'bold');

    % This is in figure coordinates, so it is always visible.
    annotationText = sprintf('Max = %.4g %s\nMin = %.4g %s', ...
                             stressMax, unitText, stressMin, unitText);
    annotation('textbox',[0.73 0.78 0.18 0.10], ...
               'String',annotationText, ...
               'FitBoxToText','on', ...
               'BackgroundColor','w', ...
               'EdgeColor','k', ...
               'FontSize',10, ...
               'FontWeight','bold');

    legend([h_def,h_und], ...
        {sprintf('Deformed stress contour (x%.2g)',scaleFactor),'Undeformed'}, ...
        'Location','best');
else
    legend([h_def,h_und], ...
        {sprintf('Deformed (x%.2g)',scaleFactor),'Undeformed'}, ...
        'Location','best');
end

%% Title
if isempty(userTitle)
    if hasStress
        titleStr = sprintf('Deformed Shape with Stress Contour (scale x%.2g)',scaleFactor);
    else
        titleStr = sprintf('Deformed vs Undeformed Shape (scale x%.2g)',scaleFactor);
    end
else
    titleStr = userTitle;
end

title(titleStr,'FontSize',13,'FontWeight','bold');
hold off;
drawnow;
end