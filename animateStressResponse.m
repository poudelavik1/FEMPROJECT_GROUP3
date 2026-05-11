%% animateStressResponse.m
% Animates stress response step by step and saves as GIF.
%
% THREE panels:
%   Left   : 3-D deformed mesh coloured by stress ratio (sigma/sigma_max)
%   Centre : vertical stress capacity bar (peak element stress %)
%   Right  : stress-displacement curve growing in real time
%
% USAGE (paste at the end of mainfile_displacement_control_FIXED.m):
%
%   animateStressResponse(XYZCoord, ELEMCon, [], ...
%       U_history, R_history, sigma_max, ...
%       'MonitorNode', loadNodes(end), ...
%       'MonitorDOF',  1, ...
%       'SaveGif',     'stress_response.gif', ...
%       'E0',          E0, 'nu', nu);

function animateStressResponse(XYZCoord, ELEMCon, ELEMENT_history, ...
                                U_history, R_history, sigma_max, varargin)

%% ── Parse inputs ─────────────────────────────────────────────────────────
p = inputParser;
addParameter(p, 'SaveGif',     'stress_response.gif');
addParameter(p, 'FrameDelay',  0.15);
addParameter(p, 'GifQuality',  128);
addParameter(p, 'MonitorNode', []);
addParameter(p, 'MonitorDOF',  1);
addParameter(p, 'ScaleFactor', []);
addParameter(p, 'E0',          30e9);
addParameter(p, 'nu',          0.18);
parse(p, varargin{:});

gifFile    = p.Results.SaveGif;
frameDelay = p.Results.FrameDelay;
gifQuality = min(max(round(p.Results.GifQuality),2),256);
monNode    = p.Results.MonitorNode;
monDOF     = p.Results.MonitorDOF;
scaleFixed = p.Results.ScaleFactor;
E0         = p.Results.E0;

nNode  = size(XYZCoord,1);
NE     = size(ELEMCon,1);
nSteps = size(U_history,2);

%% ── Monitor node (for displacement axis) ────────────────────────────────
if isempty(monNode)
    uFinal   = U_history(:,end);
    compDisp = uFinal(monDOF:3:end);
    [~, monNode] = max(abs(compDisp));
end
monGDOF = 3*(monNode-1) + monDOF;
u_curve = U_history(monGDOF,:) * 1e3;     % [mm]  — one value per step
F_curve = abs(R_history(:)') / 1e3;       % [kN]  — kept for title only

%% ── Scale factor ─────────────────────────────────────────────────────────
if isempty(scaleFixed)
    uF  = U_history(:,end);
    UX  = uF(1:3:end); UY = uF(2:3:end); UZ = uF(3:3:end);
    mxD = max(sqrt(UX.^2+UY.^2+UZ.^2));
    mxL = max(max(XYZCoord,[],1)-min(XYZCoord,[],1));
    scaleFactor = (mxD > 1e-20)*0.10*mxL/mxD + (mxD <= 1e-20);
else
    scaleFactor = scaleFixed;
end
fprintf('[animateStressResponse]  Scale = %.2f  |  Steps = %d\n',...
        scaleFactor, nSteps);

%% ── Hex8 face list ───────────────────────────────────────────────────────
hexFaces = [1 2 3 4; 5 6 7 8; 1 2 6 5; 2 3 7 6; 3 4 8 7; 4 1 5 8];
allFaces = zeros(NE*6,4,'uint32');
k = 0;
for e = 1:NE
    c = ELEMCon(e,:);
    for f = 1:6; k=k+1; allFaces(k,:)=c(hexFaces(f,:)); end
end

%% ── Axis limits ──────────────────────────────────────────────────────────
uF   = U_history(:,end);
UXf  = uF(1:3:end); UYf = uF(2:3:end); UZf = uF(3:3:end);
XYZf = XYZCoord + scaleFactor*[UXf,UYf,UZf];
allP = [XYZCoord; XYZf];
lo   = min(allP,[],1) - 0.12*max(max(allP,[],1)-min(allP,[],1));
hi   = max(allP,[],1) + 0.12*max(max(allP,[],1)-min(allP,[],1));
sp   = hi - lo;

%% ── Colormap: green → yellow → red ──────────────────────────────────────
nC   = 256;
cmap = [linspace(0.05,1.00,nC)', ...
        linspace(0.60,0.05,nC)', ...
        linspace(0.05,0.02,nC)'];

%% ════════════════════════════════════════════════════════════════
%%  PRE-COMPUTE peak stress per step  [MPa]
%%  This lets us draw the full grey preview curve on the right panel.
%% ════════════════════════════════════════════════════════════════
sig_peak_MPa = zeros(1, nSteps);
H_ref = max(XYZCoord(:,3));

for step = 1:nSteps
    if ~isempty(ELEMENT_history) && numel(ELEMENT_history) >= step
        EL = ELEMENT_history{step};
        eV = arrayfun(@(e) mean(abs(EL(e).sigma_hd)), 1:NE)';
    else
        u_s     = U_history(:,step);
        UX_s    = u_s(1:3:end); UY_s = u_s(2:3:end); UZ_s = u_s(3:3:end);
        dispMag = sqrt(UX_s.^2+UY_s.^2+UZ_s.^2);
        eps_est = dispMag / max(H_ref,1e-9);
        sig_n   = E0*eps_est ./ (1 + E0*eps_est/sigma_max);
        eV = zeros(NE,1);
        for e = 1:NE
            eV(e) = mean(sig_n(ELEMCon(e,:)));
        end
    end
    sig_peak_MPa(step) = max(eV) / 1e6;
end

sigMax_MPa = max(sig_peak_MPa);
sigLim_MPa = sigma_max / 1e6;

%% ════════════════════════════════════════════════════════════════
%%  FIGURE  — fixed axis positions, colorbar added last
%% ════════════════════════════════════════════════════════════════
fig = figure('Color','w','Name','Stress-Displacement Animation');
fig.Position = [40 40 1380 600];

ax3d  = axes('Parent',fig,'Position',[0.03 0.12 0.40 0.82]);
axBar = axes('Parent',fig,'Position',[0.50 0.12 0.09 0.72]);
axSD  = axes('Parent',fig,'Position',[0.63 0.12 0.35 0.82]);

%% ── ax3d: 3-D deformed mesh ──────────────────────────────────────────────
patch('Parent',ax3d,'Vertices',XYZCoord,'Faces',allFaces,...
      'FaceColor',[0.82 0.82 0.82],'FaceAlpha',0.15,...
      'EdgeColor',[0.65 0.65 0.65],'LineStyle','--','LineWidth',0.4);

hp = patch('Parent',ax3d,'Vertices',XYZCoord,'Faces',allFaces,...
           'FaceVertexCData',zeros(nNode,1),'FaceColor','interp',...
           'EdgeColor',[0.10 0.20 0.42],'LineWidth',0.6,'FaceAlpha',0.90);

colormap(ax3d, cmap);
clim(ax3d,[0 1]);

cb = colorbar(ax3d,'eastoutside');
cb.Label.String  = '\sigma_{sat} / \sigma_{max}';
cb.Label.FontSize = 9;
cb.Ticks      = [0 0.25 0.5 0.75 1];
cb.TickLabels = {'0%','25%','50%','75%','100%'};
cb.FontSize   = 8;
ax3d.Position = [0.03 0.12 0.40 0.82];   % restore after colorbar shift

set(ax3d,'XLim',[lo(1),hi(1)],'YLim',[lo(2),hi(2)],'ZLim',[lo(3),hi(3)]);
pbaspect(ax3d,[sp(1),sp(2),sp(3)]);
view(ax3d,3); grid(ax3d,'on');
xlabel(ax3d,'X [m]','FontSize',10);
ylabel(ax3d,'Y [m]','FontSize',10);
zlabel(ax3d,'Z [m]','FontSize',10);
set(ax3d,'FontSize',9,'Box','on');
ht3d = title(ax3d,'Step 0','FontSize',11,'FontWeight','bold');

%% ── axBar: vertical stress capacity bar ──────────────────────────────────
set(axBar,'XLim',[0 1],'YLim',[0 1],'XTick',[],...
          'YTick',0:0.25:1,...
          'YTickLabel',{'0%','25%','50%','75%','100%'},...
          'FontSize',9,'YAxisLocation','right','Box','on');
grid(axBar,'on');
ylabel(axBar,'Peak stress  /  \sigma_{max}','FontSize',9);
title(axBar,sprintf('\\sigma_{max}\n= %.0f MPa',sigLim_MPa),...
      'FontSize',9,'FontWeight','bold');

patch('Parent',axBar,'XData',[0 1 1 0],'YData',[0 0 1 1],...
      'FaceColor',[0.92 0.92 0.92],'EdgeColor','none');

hBar = patch('Parent',axBar,...
             'XData',[0 1 1 0],'YData',[0 0 1e-4 1e-4],...
             'FaceColor',[0.10 0.65 0.10],...
             'EdgeColor',[0.05 0.35 0.05],'LineWidth',1.2);

line('Parent',axBar,'XData',[0 1],'YData',[1 1],...
     'Color','r','LineStyle','--','LineWidth',1.8);
text(0.05,1.04,'\sigma_{max}','Color','r','FontSize',8.5,...
     'FontWeight','bold','Parent',axBar);

hBarTxt = text(0.5,0.5,'0%','Parent',axBar,...
               'HorizontalAlignment','center','VerticalAlignment','middle',...
               'FontSize',13,'FontWeight','bold','Color',[0.2 0.2 0.2]);

%% ── axSD: stress-displacement curve ─────────────────────────────────────
% Full preview in grey
line('Parent',axSD,'XData',u_curve,'YData',sig_peak_MPa,...
     'Color',[0.82 0.82 0.82],'LineWidth',2,'LineStyle','-');

% sigma_max capacity line
line('Parent',axSD,'XData',[min(u_curve) max(abs(u_curve))*1.1],...
     'YData',[sigLim_MPa sigLim_MPa],...
     'Color',[0.85 0 0],'LineStyle','--','LineWidth',1.5);
text(min(u_curve)+0.3, sigLim_MPa+0.5,...
     sprintf('\\sigma_{max} = %.0f MPa', sigLim_MPa),...
     'FontSize',8.5,'Color',[0.85 0 0],'FontWeight','bold','Parent',axSD);

% Peak stress marker
[pkS, pkSI] = max(sig_peak_MPa);
line('Parent',axSD,'XData',u_curve(pkSI),'YData',pkS,...
     'Marker','v','MarkerSize',9,'MarkerFaceColor',[0.85 0 0],...
     'Color',[0.85 0 0],'LineStyle','none');
text(u_curve(pkSI)+0.3, pkS+0.4,...
     sprintf('Peak = %.1f MPa', pkS),...
     'FontSize',8.5,'Color',[0.85 0 0],'FontWeight','bold','Parent',axSD);

% Growing curve
hSD = line('Parent',axSD,'XData',NaN,'YData',NaN,...
           'Color',[0.13 0.55 0.13],'LineStyle','-','Marker','o',...
           'LineWidth',2,'MarkerSize',3.5,'MarkerFaceColor',[0.13 0.55 0.13]);

% Current point
hPt = line('Parent',axSD,'XData',NaN,'YData',NaN,...
           'Marker','s','MarkerSize',11,'LineStyle','none',...
           'MarkerFaceColor',[1 0.2 0.2],'MarkerEdgeColor','w',...
           'LineWidth',1.5);

xMin = min(u_curve);
xMax = max(abs(u_curve));
yMax = max(sigLim_MPa, sigMax_MPa) * 1.18 + 1;
set(axSD,'XLim',[xMin*1.08-0.5, xMax*1.08+0.5],...
         'YLim',[-0.5, yMax],...
         'FontSize',9,'Box','on');
grid(axSD,'on');
xlabel(axSD,'Applied displacement  u  [mm]','FontSize',10);
ylabel(axSD,'Peak element stress  [MPa]','FontSize',10);
title(axSD,'Stress-Displacement Response','FontSize',11,'FontWeight','bold');

% Phase labels
midLeft  = u_curve(max(1, round(pkSI*0.40)));
if pkSI < nSteps-2
    midRight = u_curve(min(nSteps, round((pkSI+nSteps)/2)));
    text(midRight, pkS*0.50,'Softening phase','FontSize',8,...
         'Color',[0.70 0.40 0.0],'FontAngle','italic','Parent',axSD);
end
text(midLeft, pkS*0.30,'Elastic phase','FontSize',8,...
     'Color',[0.15 0.55 0.15],'FontAngle','italic','Parent',axSD);

hStepTxt = text(0.04,0.96,'','Units','normalized','FontSize',9,...
                'Color',[0.15 0.25 0.55],'FontWeight','bold','Parent',axSD);

hValTxt  = text(0.96,0.10,'','Units','normalized','FontSize',9,...
                'Color',[0.20 0.20 0.20],'HorizontalAlignment','right',...
                'Parent',axSD,...
                'BackgroundColor','w','EdgeColor',[0.80 0.80 0.80],...
                'Margin',3);

%% ── GIF init ─────────────────────────────────────────────────────────────
saveGif = ischar(gifFile) && ~isempty(gifFile);
if saveGif
    fprintf('[animateStressResponse]  Saving to %s\n', gifFile);
end

%% ════════════════════════════════════════════════════════════════
%%  ANIMATION LOOP
%% ════════════════════════════════════════════════════════════════
for step = 1:nSteps

    u_s = U_history(:,step);
    UX  = u_s(1:3:end); UY = u_s(2:3:end); UZ = u_s(3:3:end);
    XYZ_s = XYZCoord + scaleFactor*[UX,UY,UZ];

    %% Nodal stress ratio
    if ~isempty(ELEMENT_history) && numel(ELEMENT_history) >= step
        EL = ELEMENT_history{step};
        eV = arrayfun(@(e) mean(abs(EL(e).sigma_hd)), 1:NE)';
    else
        dispMag = sqrt(UX.^2+UY.^2+UZ.^2);
        eps_est = dispMag / max(H_ref,1e-9);
        sig_n   = E0*eps_est ./ (1 + E0*eps_est/sigma_max);
        eV = zeros(NE,1);
        for e = 1:NE
            eV(e) = mean(sig_n(ELEMCon(e,:)));
        end
    end

    nSum = zeros(nNode,1); nCnt = zeros(nNode,1);
    for e = 1:NE
        c = ELEMCon(e,:);
        nSum(c) = nSum(c) + eV(e);
        nCnt(c) = nCnt(c) + 1;
    end
    nCnt(nCnt==0) = 1;
    nodalRatio = min(nSum./nCnt / sigma_max, 1);
    peakRatio  = min(max(eV)/sigma_max, 1.0);
    barClr     = stressColor(peakRatio);

    %% Update 3-D patch
    set(hp,'Vertices',XYZ_s,'FaceVertexCData',nodalRatio);

    %% Update 3-D title
    set(ht3d,'String',sprintf(...
        'Step %d/%d   u = %.1f mm   \x03c3_{peak} = %.1f MPa',...
        step, nSteps, u_curve(step), sig_peak_MPa(step)));

    %% Update stress bar
    yTop = max(peakRatio, 1e-4);
    set(hBar,'XData',[0 1 1 0],'YData',[0 0 yTop yTop],...
             'FaceColor',barClr);
    set(hBarTxt,'Position',[0.5, max(yTop*0.5,0.07)],...
                'String',sprintf('%.0f%%', peakRatio*100),...
                'Color',barTextColor(peakRatio));

    %% Update stress-displacement curve
    set(hSD,'XData',u_curve(1:step),'YData',sig_peak_MPa(1:step));
    set(hPt,'XData',u_curve(step), 'YData',sig_peak_MPa(step));
    set(hStepTxt,'String',sprintf('Step %d / %d', step, nSteps));
    set(hValTxt,'String', sprintf('u = %.1f mm\n\x03c3 = %.2f MPa\n(%.0f%% of \x03c3_{max})',...
        u_curve(step), sig_peak_MPa(step), peakRatio*100));

    drawnow limitrate;

    %% Write GIF frame
    if saveGif
        frame    = getframe(fig);
        im       = frame2im(frame);
        [ind,cm] = rgb2ind(im, gifQuality);
        if step == 1
            imwrite(ind,cm,gifFile,'gif','Loopcount',inf,'DelayTime',frameDelay);
        else
            imwrite(ind,cm,gifFile,'gif','WriteMode','append','DelayTime',frameDelay);
        end
        fprintf('  Frame %2d / %d written\n', step, nSteps);
    else
        pause(frameDelay);
    end
end

%% Final hold frame
if saveGif
    frame    = getframe(fig);
    im       = frame2im(frame);
    [ind,cm] = rgb2ind(im, gifQuality);
    imwrite(ind,cm,gifFile,'gif','WriteMode','append','DelayTime',1.5);
    fprintf('[animateStressResponse]  Done → %s\n', gifFile);
end

end  % function

%% ── LOCAL helpers ────────────────────────────────────────────────────────
function c = stressColor(r)
r = max(min(r,1),0);
c = [min(2*r,1.0),  min(2*(1-r),1.0)*0.80,  0.04];
end

function c = barTextColor(r)
c = [1 1 1]*(r > 0.20) + [0.2 0.2 0.2]*(r <= 0.20);
end