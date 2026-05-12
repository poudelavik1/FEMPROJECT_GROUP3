%% MAIN_PortalFrame_Nonlinear.m
% CE 72.12 FEM Course Project 2026 – Team 3
% Portal frame under monotonic lateral load | Hex8 elements | Nonlinear material

close all; clear; clc

%% =========================
%  PATH + INPUT
%% =========================
addpath(genpath(pwd))

mshFile   = 'GMSH/frame.msh';
excelFile = 'frame_output.xlsx';

fprintf('Reading mesh from: %s\n', mshFile);
[XYZCoord, ELEMCon, NODE, ELEMENT] = Hex8_ExtractMeshToExcel(mshFile, excelFile);

nNode = size(XYZCoord,1);
NE    = size(ELEMCon,1);
nDOF  = 3*nNode;

fprintf('Nodes: %d  |  Elements: %d  |  DOFs: %d\n', nNode, NE, nDOF);
 
%% =========================
%  MATERIAL  (HD model)
%  sigma = E0*eps / (1 + E0*eps/sigma_max)
%% =========================
E0        = 30e9;    % Pa   initial Young's modulus  (30 GPa)
nu        = 0.18;    % -    Poisson's ratio
sigma_max = 35e6;    % Pa   compressive strength f'c  (30 MPa)
eps_0     = 0.002;   % -    strain at peak stress  (Hognestad)
eps_u     = 0.0035;  % -    ultimate crushing strain

%% =========================
%  BOUNDARY CONDITIONS
%% =========================
nodeTol    = 1e-3;
fixedNodes = find(abs(XYZCoord(:,3)) <= nodeTol);   % all nodes at Z = 0

fixedDOF = [];
for i = 1:numel(fixedNodes)
    n = fixedNodes(i);
    fixedDOF = [fixedDOF, 3*n-2, 3*n-1, 3*n];
end
fixedDOF = unique(fixedDOF);

fprintf('Fixed nodes : %d\n', numel(fixedNodes));
fprintf('Fixed DOFs  : %d\n', numel(fixedDOF));

%% =========================
%  LOAD VECTOR
%  Lateral load in X direction at top node(s)
%% =========================
loadNodes  = [20 95 24 176 318 179 12 98 16];
Fmax_total = -3000000;          % N  total reference load
F_total    = zeros(nDOF, 1);

F_per_node = Fmax_total / numel(loadNodes);   % -403400/9 = -44822.2 N each

for i = 1:numel(loadNodes)
    n = loadNodes(i);
    F_total(3*n - 2) = F_per_node;   % X-direction DOF of node n
end

fprintf('Load nodes   : %d\n', numel(loadNodes));
fprintf('F per node   : %.2f N  (%.2f kN)\n', F_per_node, F_per_node/1e3);
fprintf('Total F in X : %.2f N  (%.2f kN)\n', sum(F_total(1:3:end)), sum(F_total(1:3:end))/1e3);

plotmesh(excelFile, XYZCoord, fixedNodes, loadNodes)

%% =========================
%  LOAD STEPPING
%% =========================
nSteps  = 100;     %load increment per step
maxIter = 50;      % NR max iterations per step
nrTol   = 1e-3;   % NR convergence tolerance

fprintf('Load steps : %d\n', nSteps);

%% =========================
%  DISPLACEMENT STOP SETTINGS
%% =========================
dispLimit_mm = 58;                 % target displacement in mm
dispLimit    = dispLimit_mm/1000;  % convert mm to m

monitorNode = [20 95 24 176 318 179 12 98 16];
monitorDir  = 1;                   % 1 = X, 2 = Y, 3 = Z
monitorDOF  = 3*(monitorNode - 1) + monitorDir;

fprintf('Displacement stop limit : %.1f mm\n', dispLimit_mm);
fprintf('Monitor node            : %d\n', monitorNode);
fprintf('Monitor direction       : X\n');

%% =========================
%  NONLINEAR SOLVER
%% =========================
[U_history, F_history, ELEMENT, failStep, failReason] = NonlinearFEM_Solver( ...
    XYZCoord, ELEMCon, NODE, ELEMENT, ...
    E0, nu, sigma_max, ...
    F_total, fixedDOF, ...
    Fmax_total, nSteps, maxIter, nrTol);
fprintf('\nNonlinear analysis completed.\n');
fprintf('Failure info : %s\n', failReason);

%% =========================
%  POST-PROCESS
%% =========================
u_final = U_history(:, end);

UX = u_final(1:3:end);
UY = u_final(2:3:end);
UZ = u_final(3:3:end);


%% ---stress and tangent modulus summary ----------------------
sigma_hd_elem = zeros(NE,1);
E_tan_elem    = zeros(NE,1);

for eNo = 1:NE
    sigma_hd_elem(eNo) = mean(ELEMENT(eNo).sigma_hd);   % avg over 8 GPs [Pa]
    E_tan_elem(eNo)    = mean(ELEMENT(eNo).E_tan);       % avg over 8 GPs [Pa]
end

fprintf('\n--- Hognestad stress summary at final load step ---\n');
fprintf('  Max  = %.3f MPa\n', max(sigma_hd_elem)  / 1e6);
fprintf('  Mean = %.3f MPa\n', mean(sigma_hd_elem) / 1e6);
fprintf('  Min  = %.3f MPa\n', min(sigma_hd_elem)  / 1e6);
fprintf('  f''c  = %.1f MPa   (failure threshold)\n', sigma_max / 1e6);

fprintf('\n--- Tangent modulus (softening indicator) ---\n');
fprintf('  Max  = %.4f GPa\n', max(E_tan_elem)  / 1e9);
fprintf('  Mean = %.4f GPa\n', mean(E_tan_elem) / 1e9);
fprintf('  Min  = %.4f GPa  (%.1f%% of E0)\n', ...
        min(E_tan_elem)/1e9, 100*min(E_tan_elem)/E0);

%% --- Failure summary ---------------------------------------------------
if ~isempty(failStep)
    fprintf('\n*** STRUCTURE FAILED AT STEP %d ***\n', failStep);
    fprintf('%s\n', failReason);
else
    fprintf('\nStructure did not fail within %d load steps.\n', nSteps);
end

%% --- Element 1 diagnostic print ----------------------------------------
fprintf('\n--- Element 1 results ---\n');
disp('Strain at Gauss points [6 components]:');
disp(ELEMENT(1).strain);

disp('Stress at Gauss points [Pa]:');
disp(ELEMENT(1).stress);

disp('Hognestad sigma_hd at Gauss points [Pa]:');
disp(ELEMENT(1).sigma_hd);

disp('Tangent modulus E_tan at Gauss points [Pa]:');
disp(ELEMENT(1).E_tan);
%% =========================
%  DEFORMED SHAPE PLOT
%  Undeformed (grey) + deformed (blue) — no contours
%% =========================
plotDeformedShape(XYZCoord, ELEMCon, u_final);

%% =========================
%  STRESS CONTOUR PLOTS  (Hognestad model fields only)
%% =========================

% Hognestad stress — primary failure indicator
plotDeformedMesh(XYZCoord, ELEMCon, u_final, ELEMENT, ...
    'Field',  'sigma_hd', ...
    'fc',     sigma_max);
%% =========================
%  LOAD STEP ANIMATION
%% =========================
animateLoadSteps(XYZCoord, ELEMCon, U_history, F_history, ...
    'MonitorNode', loadNodes(end), ...
    'MonitorDOF',  1, ...       % 1=X  2=Y  3=Z
    'FrameDelay',  0.15);
%% Total force at loaded nodes — computed directly from F_total scaling
nStepsActual = size(U_history, 2);          % may be less than nSteps if stopped early

loadDOFs = 3*(loadNodes - 1) + 1;           % X-direction DOFs of all loaded nodes
F_ref    = abs(sum(F_total(loadDOFs)));     % total reference force magnitude [N]

F_nodes  = zeros(1, nStepsActual);
for step = 1:nStepsActual
    F_nodes(step) = F_ref * step / nSteps;  % force applied at this step [N]
end

%% Print summary
fprintf('\n--- Force at loaded nodes ---\n');
fprintf('  Loaded nodes   : %s\n', num2str(loadNodes));
fprintf('  Loaded DOFs    : %s\n', num2str(loadDOFs));
fprintf('  F per node     : %.2f kN\n', F_ref/numel(loadNodes)/1e3);
fprintf('  F total ref    : %.2f kN\n', F_ref/1e3);

