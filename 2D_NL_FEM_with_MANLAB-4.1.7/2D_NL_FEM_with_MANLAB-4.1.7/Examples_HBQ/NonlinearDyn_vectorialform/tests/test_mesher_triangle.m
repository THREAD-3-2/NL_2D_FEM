clear
close all
clc

%% Input definition
% geometric inputs (optional if mesh if directly provided)
% geometry for a triangle (length in mm)
geom_node = 1000* [1 1  0 0;              % Point 1 coord (1,0)
                   2 -1/2 sqrt(3/2) 0;    % Point 2 coord (-1/2,sqrt(3/2))
                   3 -1/2 -sqrt(3/2) 0 ;];% Point 3 coord (-1/2,-sqrt(3/2))         

geom_element = [1 1 2;  % line 1 between point 1 and 2
                2 2 3;  % line 2 between point 2 and 3
                3 3 1]; % line 3 between point 3 and 1
discretisation = [4,... % line 1 discretized with 5 elements
                  4,... % line 2 discretized with 6 elements
                  4];   % line 3 discretized with 7 elements

% % mesh inputs (optional if geometry is provided)
% mesh for a 2 elem beam
% nodes = [1 0  0 0;    % Node 1 coord (0,0)
%          2 0  1 0;    % Node 2 coord (0,1)
%          3 0  2 0 ;];%  Node 3 coord (0,2))         
% 
% connect = [1 1  2 ;  % element 1 between node 1 and 2
%            2 2 3;];  % element 3 between node 2 and 3         

% material inputs + library of materials
b = 10; h = 10;
S = b*h; I = b*h^3/12; rho = 7.8e-9; E = 210e3; eta =0.3; k=5/6;
% Boundary condition input (Dirichlet condition only)
bc_node_list{1} = struct('node', 1,'dof', [1,2,3]);

% Forcing condition input (Dirichlet condition only)
% ponctual static forces
static_ponctual_force_node_list{1} = struct('node',10,...
                                     'dof', [1 2],...
                                     'amplitude', [0 0] );
% distrubuted static forces
static_distributed_force_amplitude = [0  0];
static_distributed_force_direction = [1  2];

% ponctual periodic forces
periodic_ponctual_force_node_list{1} = struct('node',2,...
                                     'dof', [1 2],...
                                     'amplitude', [1+1i 1+2*1i],... % complex amplitude re cos + im sin 
                                     'harmonic', [1 1] );
% distrubuted periodic forces
periodic_distributed_force_amplitude = [10+2*1i 10]; % complex amplitude
periodic_distributed_force_direction = [1  2];
periodic_distributed_force_harmonic = [1  1];

%% FEM model definition
model = NL_2D_FEM; % create a new NL_2D_FEM model
% set geometry (optional if a mesh is directly provided)
model = model.set_geom(geom_node, geom_element, discretisation); % set geometry
% set mesh 
model = model.auto_mesh(); % auto mesh if geometry is given
% model = model.set_mesh(nodes, connect; % user defined mesh otherwise: needs nodes coordiantes and connectivity table
% set properties
model = model.set_prop(S, I, rho, E, eta, k); % set propoerties (must be the same for all elements)
% set boundary condition
model = model.set_boundary(bc_node_list);
% set loads (forcing)
% static loads
model = model.set_static_loads('ponctual', static_ponctual_force_node_list);
model = model.set_static_loads('distributed', static_distributed_force_amplitude, static_distributed_force_direction);
% dynamic loads
model = model.set_periodic_loads('ponctual', periodic_ponctual_force_node_list);
model = model.set_periodic_loads('distributed', periodic_distributed_force_amplitude, periodic_distributed_force_direction, periodic_distributed_force_harmonic);
% assemble mass matrix and force vector
model = model.initialise_matrices_and_vector();
% static solution
q0 = zeros(length(model.boundary.active_dof),1);
[res, fint] = model.static_residual(q0);
q0_full = zeros(3*model.mesh.number_nodes,1);
[qs_full, res] = model.solve_static_problem(q0_full);
fig = figure;
model.plot_deformed_mesh(q0_full, fig, '-k')
model.plot_deformed_mesh(qs_full, fig, '-b')
[strain,stress] = model.strains_and_stress_at_gauss_point(qs_full);

% modal analysis
[shape, freq] = model.linear_modal_analysis(qs_full);
% [shape, freq] = model.linear_modal_analysis();
fig3 = figure;
% model.plot_deformed_mesh(q0_full, fig3, '--k')
model.plot_deformed_mesh(qs_full, fig3, '--k')
model.plot_deformed_mesh(shape(:,1)*200, fig3, '-b')
model.plot_deformed_mesh(shape(:,2)*200, fig3, '-g')
model.plot_deformed_mesh(shape(:,3)*200, fig3, '-r')
% %buckling analysis
% [frq_list, lambda, N_list, q_list] = model.buckling_analysis();
% figure; plot(lambda, frq_list(1,:))
% figure; hold on; plot(lambda, q_list(4,:)); plot(lambda, q_list(5,:))
