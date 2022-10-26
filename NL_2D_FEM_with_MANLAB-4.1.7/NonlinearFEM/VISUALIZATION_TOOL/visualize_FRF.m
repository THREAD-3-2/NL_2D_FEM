function h = visualize_FRF(sys,fig,FRFs,LegendNames)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
set(0,'defaultTextInterpreter','latex');               %LaTeX pour le titre
set(0,'defaultAxesTickLabelInterpreter','latex');      %LaTeX pour les valeurs des axes
set(0,'defaultLegendInterpreter','latex');             %LaTeX pour les valeurs des légendes

Names = fieldnames(FRFs); % FRFs = struct of all forced response 'Diagram's
LineWidth = 2;
color_order = get(gca,'colororder');

obs_dof = sys.parameters.model.visu.dof;

figure(fig)
subplot(2,2,[1 3])
for jj = 1:numel(Names)
   for ii = 1:length(FRFs.(Names{jj,:}))
        h(jj) = plotbranchHBM_MPD(sys,(FRFs.(Names{jj,:}){1,ii}),obs_dof,1,'omega',LineWidth,color_order(jj,:));
        hold on
   end
   set(h(jj),'DisplayName',LegendNames(jj));
end
end

