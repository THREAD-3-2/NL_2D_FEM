function [] = point_display(sys,Uf)
% function [] = point_display(sys,Uf)
% User-defined function to display solutions along the branch.
%
% See the repertory MANLAB-4.x/SRC/Display for more informations.

fig = figure(11); fig.Name='Mesh Animation';
clf; hold on;
% Plot the variable 1,2 over one period.
%plotperiodHBM(sys,Uf,[1 2]);

%figure(12);clf; hold on;
% Plot the variable 1 with respect to variable 2 
% over one period. [phase diagram] 
%plotphasediagHBM(sys,Uf,[1 2]);

%%% Plot the amplitude of the sine and the cosine coefficients of the
%%% first and the second variable.
% plotbarsincosHBM(sys,Uf,[1 2]);

%%% Plot the amplitude of the harmonics of the first and the second
%%% variable.
% plotbarHBM(sys,Uf,[1 2])

%%% Compute the evolution in the time domain of the first and the second
%%% variable over 10 periods (the time is dimensionless) with 1e4 points.
% nb_period = 1;
% nb_samples = 1e4;
%time = linspace(0,nb_period,nb_samples);
%Utime = calcperiodHBM(sys,Uf,[1:sys.nz/2],time);
%plot(time,Utime(:,sys.parameters.dof_info.obs_dof)) %% transverse displacement at cantilever tip

% fig_number = 22;
nb_samples_animate = 16;
animate_solution(sys,Uf,fig,nb_samples_animate)
%%% Compute and plot the time-dependent displacement of the beam according to Utime
% Defining the degrees of freedom and parameters
% dof = sys.parameters.dof_info.dof;
% bc_dof = sys.parameters.dof_info.bc.prescribed_dof;
% active_dof = sys.parameters.dof_info.active;
% 
% number_nodes = sys.parameters.mesh.number_nodes;
% number_elements = sys.parameters.mesh.number_elements;
% Le = sys.parameters.beam_params.L_beam/number_elements;

end
