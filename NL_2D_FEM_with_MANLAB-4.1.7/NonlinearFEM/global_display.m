function [] = global_display(sys,Section)
% function [] = global_display(sys,Section)
% User-defined function to display the solution-branch.
%
% Section is a 'CheckPoint class' object 
% containing all the information about the series:
% discrete representation, bifurcation analysis,...
%
% See the repertory MANLAB-4.x/SRC/Display for more functions.
%%% Plot the norm of the first harmonic of the two first variables with
%%% respect to lambda
% plotbranchHBM(sys,Section,[sys.nz/2-1 sys.nz/2-1 sys.nz/2-1],[1 2 3],'omega');
% plot first harmonic Amplitude vs frequency
dispdof=sys.parameters.model.visu.dof;
for i=1:length(dispdof)
   figureFRF(i) = figure(i+200); 
   figureFRF(i).Name = ['Harmonic Amplitude vs Frequency, DOF ' num2str(dispdof(i))];  hold on;
   hand{i} = plotbranchHBM(sys,Section,dispdof(i),1,'omega');
% For rescalling grapha
%    xt = get(gca, 'XTick');                                 % 'XTick' Values
%    set(gca, 'XTick', xt, 'XTickLabel', xt/22.41)
%    yt = get(gca, 'YTick');                                 % 'XTick' Values
%    set(gca, 'YTick', yt, 'YTickLabel', yt*1000)
%    grid on
end


end
