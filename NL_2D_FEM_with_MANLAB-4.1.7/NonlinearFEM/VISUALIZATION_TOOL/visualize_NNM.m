function h = visualize_NNM(sys,fig,Diagram_NNM)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
set(0,'defaultTextInterpreter','latex');               %LaTeX pour le titre
set(0,'defaultAxesTickLabelInterpreter','latex');      %LaTeX pour les valeurs des axes
set(0,'defaultLegendInterpreter','latex');             %LaTeX pour les valeurs des légendes

LineWidth = 2;
color_order = get(gca,'colororder');

obs_dof = sys.parameters.model.visu.dof;

figure(fig)
subplot(2,2,[1 3])
hold on
for ii = 1:length(Diagram_NNM(1,:))
    h = plotbranchHBM_MPD(sys,Diagram_NNM(1,ii),obs_dof,1,'omega',LineWidth,color_order(4,:));
    hold on
end
set(h,'DisplayName','NNM');
end