function plotmesh(filename, XYZCoord, fixedNodes, loadNodes)
%PLOTMESH  Plot the portal frame mesh with fixed supports and load point.
%
%  Usage:
%    plotmesh(excelFile)
%    plotmesh(excelFile, XYZCoord, fixedNodes)
%    plotmesh(excelFile, XYZCoord, fixedNodes, loadNodes)

%% =========================================================
%  READ EXCEL
%% =========================================================
nodes  = readmatrix(filename, 'Sheet', 'Nodes');
elems  = readmatrix(filename, 'Sheet', 'Elements');

GmshID = nodes(:,2);
X      = nodes(:,3);
Y      = nodes(:,4);
Z      = nodes(:,5);

idMap  = containers.Map(GmshID, 1:length(GmshID));
memcon = elems(:,3:10);

faces = [1 2 3 4;
         5 6 7 8;
         1 2 6 5;
         2 3 7 6;
         3 4 8 7;
         4 1 5 8];

%% =========================================================
%  FIGURE SETUP
%% =========================================================
figure('Color','w','Position',[100 100 1000 700]);
hold on; axis equal; grid on; view(3);
xlabel('X [m]','FontSize',11);
ylabel('Y [m]','FontSize',11);
zlabel('Z [m]','FontSize',11);
title('Portal Frame Mesh — Boundary Conditions & Load', ...
      'FontSize',13,'FontWeight','bold');

%% =========================================================
%  DRAW MESH ELEMENTS
%% =========================================================
nElem = size(memcon,1);
for e = 1:nElem
    elemNodes = memcon(e,:);
    idx = zeros(1,8);
    for k = 1:8
        idx(k) = idMap(elemNodes(k));
    end
    vertices = [X(idx), Y(idx), Z(idx)];
    patch('Vertices', vertices, ...
          'Faces',    faces, ...
          'FaceColor','cyan', ...
          'FaceAlpha', 0.10, ...
          'EdgeColor','k', ...
          'LineWidth', 0.8);
end

%% All nodes — small grey dots
plot3(X, Y, Z, 'k.', 'MarkerSize', 4);

%% Node number labels
for i = 1:length(GmshID)
    text(X(i), Y(i), Z(i), sprintf('%d', GmshID(i)), ...
        'FontSize',           7, ...
        'Color',              [0.15 0.15 0.70], ...
        'HorizontalAlignment','left', ...
        'VerticalAlignment',  'bottom', ...
        'Clipping',           'on');
end

%% =========================================================
%  FIXED NODES — triangle markers + per-column base box
%% =========================================================
if nargin >= 3 && ~isempty(fixedNodes) && ~isempty(XYZCoord)

    Xf = XYZCoord(fixedNodes, 1);
    Yf = XYZCoord(fixedNodes, 2);
    Zf = XYZCoord(fixedNodes, 3);

    %% Triangle markers — LineStyle none prevents the connecting line
    plot3(Xf, Yf, Zf, 'r^', ...
          'LineStyle',       'none', ...   % <-- this was the bug
          'MarkerSize',       8, ...
          'LineWidth',        1.5, ...
          'MarkerFaceColor', [1 0.2 0.2]);

    %% Hatch ticks below each node
    hatchLen = 0.04 * (max(Z) - min(Z));
    for i = 1:numel(fixedNodes)
        plot3([Xf(i) Xf(i)], [Yf(i) Yf(i)], [Zf(i) Zf(i)-hatchLen], ...
              'r-', 'LineWidth', 1.2, 'HandleVisibility','off');
    end

    %% Per-column base rectangle
    %  Cluster fixed nodes by X so each column gets its OWN small box.
    %  Using a tolerance of 10% of frame width to separate clusters.
    frameWidth = max(X) - min(X);
    xTol       = 0.10 * frameWidth;

    xClusters  = uniqueClusters(Xf, xTol);   % column X centres

    for c = 1:numel(xClusters)
        inCluster = abs(Xf - xClusters(c)) <= xTol;
        xc = Xf(inCluster);  yc = Yf(inCluster);  zc = Zf(inCluster);

        pad = 0.03;
        bx  = [min(xc)-pad, max(xc)+pad, max(xc)+pad, min(xc)-pad, min(xc)-pad];
        bz  = [min(zc)-pad, min(zc)-pad, max(zc)+pad, max(zc)+pad, min(zc)-pad];
        by  = min(yc) * ones(1,5);

        plot3(bx, by, bz, 'r-', 'LineWidth', 2.0, 'HandleVisibility','off');

        %% FIXED BASE label centered under each column
        text(mean(xc), min(yc), min(zc) - 2*hatchLen, ...
            'FIXED', ...
            'FontSize',           8, ...
            'FontWeight',         'bold', ...
            'Color',              [0.8 0 0], ...
            'HorizontalAlignment','center');
    end
end

%% =========================================================
%  LOAD NODES — arrow + label
%% =========================================================
if nargin >= 4 && ~isempty(loadNodes) && ~isempty(XYZCoord)

    arrowLen = 0.15 * (max(X) - min(X));

    for k = 1:numel(loadNodes)
        n  = loadNodes(k);
        Xl = XYZCoord(n,1);
        Yl = XYZCoord(n,2);
        Zl = XYZCoord(n,3);

        %% Star marker at load node
        plot3(Xl, Yl, Zl, 'p', ...
              'MarkerSize',      16, ...
              'LineWidth',        2, ...
              'MarkerFaceColor', [0 0.8 0], ...
              'MarkerEdgeColor', [0 0.5 0]);

        %% Arrow in -X direction
        quiver3(Xl + arrowLen, Yl, Zl, ...
                -arrowLen, 0, 0, ...
                0, ...
                'Color',      [0 0.6 0], ...
                'LineWidth',   2.5, ...
                'MaxHeadSize', 0.8, ...
                'HandleVisibility','off');

        %% Label
        text(Xl + arrowLen + 0.02, Yl, Zl, ...
            sprintf(' F\n Node %d\n (%.2f,%.2f,%.2f)', n, Xl, Yl, Zl), ...
            'FontSize',           9, ...
            'FontWeight',         'bold', ...
            'Color',              [0 0.55 0], ...
            'VerticalAlignment',  'middle', ...
            'BackgroundColor',    [0.9 1 0.9], ...
            'EdgeColor',          [0 0.6 0], ...
            'Margin',             3);
    end
end

%% =========================================================
%  LEGEND
%% =========================================================
hMesh = plot3(NaN,NaN,NaN,'c-','LineWidth',2,'DisplayName','Mesh edges');
legendEntries = {hMesh};

if nargin >= 3 && ~isempty(fixedNodes)
    hFix = plot3(NaN,NaN,NaN,'r^', ...
                 'MarkerFaceColor','r','MarkerSize',8, ...
                 'LineStyle','none','DisplayName','Fixed nodes');
    legendEntries{end+1} = hFix;
end

if nargin >= 4 && ~isempty(loadNodes)
    hLoad = plot3(NaN,NaN,NaN,'p', ...
                  'MarkerFaceColor',[0 0.8 0], ...
                  'MarkerEdgeColor',[0 0.5 0], ...
                  'MarkerSize',12,'DisplayName','Load node');
    legendEntries{end+1} = hLoad;
end

legend([legendEntries{:}],'Location','best','FontSize',10);
hold off;
drawnow;
end

%% =========================================================
%  LOCAL: find cluster centres given a tolerance
%% =========================================================
function centres = uniqueClusters(vals, tol)
    sorted  = sort(vals);
    centres = sorted(1);
    for i = 2:numel(sorted)
        if abs(sorted(i) - centres(end)) > tol
            centres(end+1) = sorted(i); %#ok<AGROW>
        end
    end
end