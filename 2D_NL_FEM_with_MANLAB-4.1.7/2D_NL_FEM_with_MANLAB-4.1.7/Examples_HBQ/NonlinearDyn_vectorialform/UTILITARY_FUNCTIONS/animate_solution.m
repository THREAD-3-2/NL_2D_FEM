function animate_solution(sys,U,figure_number,nb_samples)
% Plot the time evolution of the mesh 

set(0,'defaultTextInterpreter','latex');               %LaTeX pour le titre
set(0,'defaultAxesTickLabelInterpreter','latex');      %LaTeX pour les valeurs des axes
set(0,'defaultLegendInterpreter','latex');             %LaTeX pour les valeurs des légendes

% Defining the degrees of freedom and parameters
model = sys.parameters.model;
active_dof = model.boundary.active_dof;
number_nodes = model.mesh.number_nodes;

nb_period = 1;
% nb_samples = 1.6^6;
time = linspace(0,nb_period,nb_samples+1);
Utime = calcperiodHBM(sys,U,[1:sys.nz/2],time);


fig = figure(figure_number);

% Initialization
q_global = zeros(3*number_nodes,size(Utime,1));
% find mesh limits
maxXYmesh = max(model.mesh.nodes(:,2:3));
minXYmesh = min(model.mesh.nodes(:,2:3));
% find displacement limits
q_global(active_dof,:) = Utime';
maxXdisp = max(max(abs(q_global(1:3:end,:))));
maxYdisp = max(max(abs(q_global(2:3:end,:))));
% define axis limits
axis([(minXYmesh(1)-maxXdisp)*1.1, (maxXYmesh(1)+maxXdisp)*1.1, (minXYmesh(2)-maxYdisp)*1.1, (maxXYmesh(2)+maxYdisp)*1.1])
for j = 1:length(time)
    q_global(active_dof) = Utime(j,:);
%cla
    model.plot_deformed_mesh(q_global, fig,'-b') % plot deformed mesh at time tj
    drawnow % forces to draw all figure now
    pause(0.05)
end
model.plot_deformed_mesh(0*q_global, fig,'--k') % plot undeformed mesh 
end



