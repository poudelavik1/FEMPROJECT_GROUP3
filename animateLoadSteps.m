function animateLoadSteps(XYZCoord, ELEMCon, U_history, F_history, varargin)
%ANIMATELOADSTEPS  Animate deformed mesh through each load step.
%
%  Shows the portal frame deforming as load increases from 0 to failure.
%  Left panel : 3-D deformed mesh coloured by displacement magnitude.
%  Right panel: growing load-displacement curve with current point.
%
%  Usage
%    animateLoadSteps(XYZCoord, ELEMCon, U_history, F_history)
%    animateLoadSteps(..., 'MonitorNode', 318, 'MonitorDOF', 1)
%    animateLoadSteps(..., 'SaveGif', 'portal_frame.gif')
%
%  Options (name-value)
%    'MonitorNode'  — node number to track on load-disp curve (default: last loaded node)
%    'MonitorDOF'   — local DOF of monitor node: 1=X 2=Y 3=Z  (default: 1)
%    'ScaleFactor'  — deformation magnification ([] = auto per step)
%    'FrameDelay'   — pause between frames in seconds (default: 0.08)
%    'SaveGif'      — filename to save animated GIF, e.g. 'anim.gif'
%                     leave empty '' to not save (default: '')
%    'Fmax'         — reference load for % label [N] (default: max F_history)

%% ========================================================================
%  Parse inputs
%% ========================================================================
p = inputParser;
addParameter(p, 'MonitorNode',  []);
addParameter(p, 'MonitorDOF',   1);
addParameter(p, 'ScaleFactor',  []);
addParameter(p, 'FrameDelay',   0.08);
addParameter(p, 'SaveGif',      '');
addParameter(p, 'Fmax',         []);
parse(p, varargin{:});

monNode    = p.Results.MonitorNode;
monDOF     = p.Results.MonitorDOF;
scaleFixed = p.Results.ScaleFactor;
frameDelay = p.Results.FrameDelay;
gifFile    = p.Results.SaveGif;
Fmax       = p.Results.Fmax;

nNode  = size(XYZCoord, 1);
NE     = size(ELEMCon,  1);
nSteps = size(U_history, 2);

if isempty(Fmax)
    Fmax = max(abs(F_history));
end

%% Monitor DOF for load-disp curve
if isempty(monNode)
    % Default: use the node with highest final displacement in monitored direction
    uFinal = U_history(:, end);
    compDisp = uFinal(monDOF:3:end);
    [~, monNode] = max(abs(compDisp));
end
monGlobalDOF = 3*(monNode-1) + monDOF;

%% ========================================================================
%  Hex8 face list — built once, reused every frame
%% ========================================================================
hexFaces = [1 2 3 4;
            5 6 7 8;
            1 2 6 5;
            2 3 7 6;
            3 4 8 7;
            4 1 5 8];

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
%  Global scale factor — computed from final step so it stays constant
%% ========================================================================
if isempty(scaleFixed)
    uFinal = U_history(:,end);
    UX = uFinal(1:3:end); UY = uFinal(2:3:end); UZ = uFinal(3:3:end);
    maxDisp = max(sqrt(UX.^2 + UY.^2 + UZ.^2));
    maxDim  = max(max(XYZCoord,[],1) - min(XYZCoord,[],1));
    if maxDisp > 1e-20
        scaleFactor = 0.10 * maxDim / maxDisp;
    else
        scaleFactor = 1;
    end
else
    scaleFactor = scaleFixed;
end
fprintf('[animateLoadSteps]  Scale factor = %.2f  |  Steps = %d\n', scaleFactor, nSteps);

%% Global colour limits — displacement magnitude at final step [mm]
uFinal = U_history(:,end);
UX = uFinal(1:3:end); UY = uFinal(2:3:end); UZ = uFinal(3:3:end);
dispMagFinal = sqrt(UX.^2 + UY.^2 + UZ.^2) * 1e3;
climMax = max(dispMagFinal);
if climMax < 1e-10; climMax = 1; end

%% Global axis limits — from final deformed + undeformed
XYZ_final = XYZCoord + scaleFactor * [UX, UY, UZ];
allPts = [XYZCoord; XYZ_final];
lo = min(allPts,[],1);  hi = max(allPts,[],1);
pad = 0.10 * max(hi - lo);
lo = lo - pad;  hi = hi + pad;

%% Load-disp curve data
u_curve = U_history(monGlobalDOF, :) * 1e3;   % [mm]
F_curve = abs(F_history) / 1e3;               % [kN]

%% ========================================================================
%  Build figure layout
%% ========================================================================
fig = figure('Color','w','Name','Load Step Animation');
fig.Position = [60 60 1200 620];

%% Left: 3-D mesh panel
ax3d = axes('Parent', fig, 'Position', [0.03 0.08 0.58 0.86]);

%% Right: load-displacement panel
axLD = axes('Parent', fig, 'Position', [0.68 0.12 0.29 0.75]);

%% ========================================================================
%  Draw initial (undeformed) state on ax3d
%% ========================================================================
axes(ax3d);

% Undeformed ghost
patch('Vertices',  XYZCoord, ...
      'Faces',     allFaces, ...
      'FaceColor', [0.80 0.80 0.80], ...
      'FaceAlpha', 0.20, ...
      'EdgeColor', [0.60 0.60 0.60], ...
      'LineStyle', '--', ...
      'LineWidth', 0.5);
hold on;

% Deformed patch — will update vertices each frame
hp = patch('Vertices',        XYZCoord, ...     % start = undeformed
           'Faces',           allFaces, ...
           'FaceVertexCData', zeros(nNode,1), ...
           'FaceColor',       'interp', ...
           'EdgeColor',       [0.10 0.22 0.45], ...
           'LineWidth',       0.7, ...
           'FaceAlpha',       0.85);

colormap(ax3d, jet(128));
cb = colorbar(ax3d, 'southoutside');
cb.Label.String = 'Displacement Magnitude [mm]';
cb.FontSize = 9;
clim(ax3d, [0, climMax]);

span3d = hi - lo;
xlim(ax3d, [lo(1), hi(1)]);
ylim(ax3d, [lo(2), hi(2)]);
zlim(ax3d, [lo(3), hi(3)]);
pbaspect(ax3d, [span3d(1), span3d(2), span3d(3)]);
grid(ax3d, 'on');
view(ax3d, 3);
xlabel(ax3d, 'X [m]','FontSize',10);
ylabel(ax3d, 'Y [m]','FontSize',10);
zlabel(ax3d, 'Z [m]','FontSize',10);
set(ax3d,'FontSize',9,'Box','on');

htitle = title(ax3d, 'Step 0 / 0  |  Load = 0.0 kN  (0%)', ...
               'FontSize',11,'FontWeight','bold');

%% Load arrow annotation (updated each step)
% We'll use a text label at the load application point instead
hLoadTxt = text(ax3d, hi(1)-pad*0.5, (lo(2)+hi(2))/2, hi(3)-pad*0.5, ...
    '', 'FontSize', 9, 'Color', [0.8 0 0], 'FontWeight', 'bold', ...
    'HorizontalAlignment', 'center');

%% ========================================================================
%  Load-displacement axes — draw full curve in grey first
%% ========================================================================
axes(axLD);
plot(axLD, u_curve, F_curve, '-', 'Color', [0.80 0.80 0.80], 'LineWidth', 1.5);
hold(axLD, 'on');

% Growing curve — plotted in blue, updated each step
hCurve = plot(axLD, NaN, NaN, 'b-o', 'LineWidth', 1.8, ...
              'MarkerSize', 4, 'MarkerFaceColor', 'b');

% Current point marker
hPt = plot(axLD, NaN, NaN, 'rs', 'MarkerSize', 10, ...
           'LineWidth', 2, 'MarkerFaceColor', [1 0.3 0.3]);

xlabel(axLD, 'Displacement [mm]','FontSize',10);
ylabel(axLD, 'Applied Load [kN]','FontSize',10);
title(axLD, 'Load-Displacement Curve','FontSize',11,'FontWeight','bold');
grid(axLD, 'on');
xlim(axLD, [min(u_curve)*1.1 - 0.01, max(u_curve)*1.1 + 0.01]);
ylim(axLD, [0, max(F_curve)*1.15 + 0.01]);
set(axLD,'FontSize',9,'Box','on');

%% Step counter text on load-disp plot
hStepTxt = text(axLD, 0.05, 0.92, '', 'Units','normalized', ...
    'FontSize', 9, 'Color', [0.2 0.2 0.6], 'FontWeight','bold');

%% ========================================================================
%  GIF initialisation
%% ========================================================================
saveGif = ~isempty(gifFile);
if saveGif
    fprintf('Saving GIF to: %s\n', gifFile);
end

%% ========================================================================
%  Animation loop
%% ========================================================================
for step = 1:nSteps

    u_step = U_history(:, step);
    UX = u_step(1:3:end);
    UY = u_step(2:3:end);
    UZ = u_step(3:3:end);

    %% Deformed coordinates
    XYZ_step = XYZCoord + scaleFactor * [UX, UY, UZ];

    %% Displacement magnitude per node [mm]
    dispMag = sqrt(UX.^2 + UY.^2 + UZ.^2) * 1e3;

    %% Update 3-D patch
    set(hp, 'Vertices',        XYZ_step, ...
            'FaceVertexCData', dispMag);

    %% Update title
    loadPct = 100 * abs(F_history(step)) / Fmax;
    set(htitle, 'String', sprintf('Step %d / %d  |  Load = %.1f kN  (%.0f%%)', ...
        step, nSteps, F_curve(step), loadPct));

    %% Update load-disp curve
    set(hCurve, 'XData', u_curve(1:step), 'YData', F_curve(1:step));
    set(hPt,    'XData', u_curve(step),   'YData', F_curve(step));
    set(hStepTxt,'String', sprintf('Step %d / %d', step, nSteps));

    drawnow;

    %% Save frame to GIF
    if saveGif
        frame  = getframe(fig);
        im     = frame2im(frame);
        [ind, cm] = rgb2ind(im, 128);
        if step == 1
            imwrite(ind, cm, gifFile, 'gif', ...
                'Loopcount', inf, 'DelayTime', frameDelay);
        else
            imwrite(ind, cm, gifFile, 'gif', ...
                'WriteMode', 'append', 'DelayTime', frameDelay);
        end
    else
        pause(frameDelay);
    end

end

%% Final frame — hold for 1 s
pause(1.0);
fprintf('[animateLoadSteps]  Done.\n');

end