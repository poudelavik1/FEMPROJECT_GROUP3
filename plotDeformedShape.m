function plotDeformedShape(XYZCoord, ELEMCon, u_final, varargin)
%PLOTDEFORMEDSHAPE  Overlay deformed and undeformed mesh — no contours.
%
%  Usage
%    plotDeformedShape(XYZCoord, ELEMCon, u_final)
%    plotDeformedShape(..., 'ScaleFactor', 200)
%
%  Options
%    'ScaleFactor'  — deformation magnification ([] = auto)
%    'Padding'      — fraction of structure size added as margin (default 0.08)
%    'Title'        — custom title string

%% Parse inputs
p = inputParser;
addParameter(p, 'ScaleFactor', []);
addParameter(p, 'Padding',     0.08);
addParameter(p, 'Title',       []);
parse(p, varargin{:});
scaleFactor = p.Results.ScaleFactor;
padding     = p.Results.Padding;
userTitle   = p.Results.Title;

%% Displacements
UX = u_final(1:3:end);
UY = u_final(2:3:end);
UZ = u_final(3:3:end);

%% Auto scale factor — deformation = 10% of bounding box
if isempty(scaleFactor)
    maxDisp = max(sqrt(UX.^2 + UY.^2 + UZ.^2));
    maxDim  = max(max(XYZCoord,[],1) - min(XYZCoord,[],1));
    if maxDisp > 1e-20
        scaleFactor = 0.10 * maxDim / maxDisp;
    else
        scaleFactor = 1;
    end
    fprintf('[plotDeformedShape]  Auto scale factor = %.2f\n', scaleFactor);
end

%% Deformed coordinates
XYZ_def = XYZCoord + scaleFactor * [UX, UY, UZ];

%% Tight axis limits — union of undeformed and deformed bounds + padding
allPts  = [XYZCoord; XYZ_def];
lo      = min(allPts, [], 1);
hi      = max(allPts, [], 1);
span    = hi - lo;
pad     = padding * max(span);   % uniform padding in metres
lo      = lo - pad;
hi      = hi + pad;

%% Hex8 face definitions (local node indices)
hexFaces = [1 2 3 4;   % bottom
            5 6 7 8;   % top
            1 2 6 5;   % front
            2 3 7 6;   % right
            3 4 8 7;   % back
            4 1 5 8];  % left

NE = size(ELEMCon, 1);
allFaces = zeros(NE*6, 4, 'uint32');
fIdx = 0;
for eNo = 1:NE
    conn = ELEMCon(eNo,:);
    for f = 1:6
        fIdx = fIdx + 1;
        allFaces(fIdx,:) = conn(hexFaces(f,:));
    end
end

%% Figure
figure('Color','w','Name','Deformed vs Undeformed');
set(gcf,'Position',[100 100 900 650]);

%% Undeformed — light grey fill, dashed edges
h_und = patch('Vertices',  XYZCoord, ...
              'Faces',     allFaces, ...
              'FaceColor', [0.80 0.80 0.80], ...
              'FaceAlpha', 0.25, ...
              'EdgeColor', [0.55 0.55 0.55], ...
              'LineStyle', '--', ...
              'LineWidth', 0.7);
hold on;

%% Deformed — solid blue fill, dark edges
h_def = patch('Vertices',  XYZ_def, ...
              'Faces',     allFaces, ...
              'FaceColor', [0.18 0.44 0.78], ...
              'FaceAlpha', 0.70, ...
              'EdgeColor', [0.05 0.20 0.45], ...
              'LineStyle', '-', ...
              'LineWidth', 0.9);
hold off;

%% Legend
legend([h_def, h_und], ...
    {sprintf('Deformed (\\times%.0f)', scaleFactor), 'Undeformed'}, ...
    'Location', 'best', 'FontSize', 10);

%% Axes — pbaspect keeps geometry correct without inflating empty space
span = hi - lo;
view(3); grid on;
xlim([lo(1), hi(1)]);
ylim([lo(2), hi(2)]);
zlim([lo(3), hi(3)]);
pbaspect([span(1), span(2), span(3)]);
xlabel('X [m]','FontSize',11);
ylabel('Y [m]','FontSize',11);
zlabel('Z [m]','FontSize',11);
set(gca,'FontSize',10,'Box','on');

if isempty(userTitle)
    titleStr = sprintf('Deformed vs Undeformed  (scale \\times%.0f)', scaleFactor);
else
    titleStr = userTitle;
end
title(titleStr,'FontSize',13,'FontWeight','bold');

drawnow;
end