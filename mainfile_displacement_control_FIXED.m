%% MAIN_PortalFrame_DisplacementControl_FIXED.m
% Displacement control version: prescribe u = 0.04h and compute reaction.
% Use this instead of force-control mainfile.m.

close all; clear; clc
addpath(genpath(pwd))

%% Mesh input
mshFile   = 'GMSH/frame.msh';
excelFile = 'frame_output.xlsx';

fprintf('Reading mesh from: %s\n', mshFile);
[XYZCoord, ELEMCon, NODE, ELEMENT] = Hex8_ExtractMeshToExcel(mshFile, excelFile);

nNode = size(XYZCoord,1);
nDOF  = 3*nNode;
fprintf('Nodes: %d | Elements: %d | DOFs: %d\n', nNode, size(ELEMCon,1), nDOF);
%% Material: Saturation model
E0        = 30e9;    % Pa
nu        = 0.18;
sigma_max = 35e6;   % Pa

%% Fixed supports: all nodes at base Z = 0
nodeTol    = 1e-6;
fixedNodes = find(abs(XYZCoord(:,3)) <= nodeTol);

fixedDOF = [];
for i = 1:numel(fixedNodes)
    n = fixedNodes(i);
    fixedDOF = [fixedDOF, 3*n-2, 3*n-1, 3*n]; 
end
fixedDOF = unique(fixedDOF);

fprintf('Fixed nodes = %d\n', numel(fixedNodes));


%% Displacement control
% Keep the same loaded node you used in the original force-control code.
loadNodes  = [318];

%DOF 1 = X
%DOF 2 = Y
%DOF 3 = Z
%X DOF = 3*n - 2
%Y DOF = 3*n - 1
%Z DOF = 3*n

controlDOF = 3*loadNodes - 2;   % X-direction DOF

h = 1.45;
u_target = -0.04*h;            % -0.058 m = -58 mm

nSteps  = 30;
maxIter = 50;
nrTol   = 1e-3;

fprintf('Target displacement = %.6f m = %.2f mm\n', u_target, u_target*1000);
plotmesh(excelFile, XYZCoord, fixedNodes, loadNodes);
%% Solve
[U_history, R_history, disp_history, ELEMENT, failStep, failReason] = ...
    NonlinearFEM_Solver_DisplacementControl_FIXED( ...
    XYZCoord, ELEMCon, NODE, ELEMENT, ...
    E0, nu, sigma_max, ...
    fixedDOF, controlDOF, u_target, ...
    nSteps, maxIter, nrTol);

fprintf('\nStatus: %s\n', failReason);
u_final = U_history(:,end);

%% =========================
% MAX/MIN STRESS + REACTIONS
%% =========================

% Element stress value from sigma_hd
stressValue = zeros(size(ELEMCon,1),1);

for eNo = 1:size(ELEMCon,1)
    stressValue(eNo) = max(ELEMENT(eNo).sigma_hd);   % Pa
end

maxStress = max(stressValue);
minStress = min(stressValue);

fprintf('\n===== STRESS RESULTS =====\n');
fprintf('Maximum stress = %.4f MPa\n', maxStress/1e6);
fprintf('Minimum stress = %.4f MPa\n', minStress/1e6);

%% Reactions at supports
% Need final internal force vector
F_int_final = localAssembleInternalForce(XYZCoord, ELEMCon, ELEMENT, nNode);

Rx = sum(F_int_final(fixedDOF(1:3:end)));
Ry = sum(F_int_final(fixedDOF(2:3:end)));
Rz = sum(F_int_final(fixedDOF(3:3:end)));

fprintf('\n===== SUPPORT REACTIONS =====\n');
fprintf('Rx = %.4f kN\n', Rx/1000);
fprintf('Ry = %.4f kN\n', Ry/1000);
fprintf('Rz = %.4f kN\n', Rz/1000);

%% =========================
% DISPLACEMENT RESULTS
%% =========================

UX = u_final(1:3:end);
UY = u_final(2:3:end);
UZ = u_final(3:3:end);

U_total = sqrt(UX.^2 + UY.^2 + UZ.^2);

maxUX = max(abs(UX));
maxUY = max(abs(UY));
maxUZ = max(abs(UZ));
maxUTotal = max(U_total);

fprintf('\n===== DISPLACEMENT RESULTS =====\n');
fprintf('Maximum X displacement = %.4f mm\n', maxUX*1000);
fprintf('Maximum Y displacement = %.4f mm\n', maxUY*1000);
fprintf('Maximum Z displacement = %.4f mm\n', maxUZ*1000);
fprintf('Maximum total displacement = %.4f mm\n', maxUTotal*1000);


%% =========================
% CONTROLLED DOF RESULTS
%% =========================

fprintf('\n===== CONTROLLED DOF RESULTS =====\n');
fprintf('Controlled nodes: ');
fprintf('%d ', loadNodes);
fprintf('\n');

fprintf('Controlled direction: X\n');
fprintf('Target displacement = %.4f mm\n', u_target*1000);
fprintf('Final controlled displacement = %.4f mm\n', disp_history(end)*1000);
fprintf('Force at controlled DOFs = %.4f kN\n', R_history(end)/1000);

%% =========================
%  DEFORMED SHAPE PLOT
%  Undeformed (grey) + deformed (blue) — no contours
%% =========================
plotDeformedMesh(XYZCoord, ELEMCon, u_final, ELEMENT, ...
    'Field',  'sigma_hd', ...
    'fc',     sigma_max);

%% =========================
%  STRESS CONTOUR PLOTS  
%% =========================
plotDeformedShape(XYZCoord, ELEMCon, u_final);

%% =========================
%  LOAD STEP ANIMATION
%% =========================
animateLoadSteps(XYZCoord, ELEMCon, U_history, abs(R_history), ...
    'MonitorNode', loadNodes(end), ...
    'MonitorDOF',  1, ...
    'FrameDelay',  0.2, ...
    'SaveGif',     'portal_frame_animation.gif', ...
    'GifQuality',  256);