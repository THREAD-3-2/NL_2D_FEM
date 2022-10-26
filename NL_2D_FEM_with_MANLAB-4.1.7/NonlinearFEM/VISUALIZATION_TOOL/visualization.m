function visualization(matFileName)
%UNTITLED Summary of this function goes here
% The order of the .mat file MUST be (the name in ' ' is what is exported
% from MANLAB window):
% FRFs 'Diagram' (1:n)
% NNM 'Diagram' (n+1)
% Legend names (saved in a vector) (n+2)
% High amplitude 'Point' (n+3)
% Low amplitude 'Point' (n+4)
set(0,'defaultTextInterpreter','latex');               %LaTeX pour le titre
set(0,'defaultAxesTickLabelInterpreter','latex');      %LaTeX pour les valeurs des axes
set(0,'defaultLegendInterpreter','latex');             %LaTeX pour les valeurs des légendes

% Diagrams = load(matFileName);
load(matFileName)

fig = figure(24);
%% Forced response curves and backbone
h1 = visualize_FRF(sys,fig,FRFs,LegendNames); % returns h1 and h2 for the legend
h2 = visualize_NNM(sys,fig,Diagram_NNM);

axis square
xlabel('$\Omega / \omega_N$')
xlim([0.75,1.1])
%% Deformed shape (low amplitude & high amplitude)
% visualize_point(sys,fig,Defs)

h = [h1 h2];
legend(h,'Location','NorthWest')
end