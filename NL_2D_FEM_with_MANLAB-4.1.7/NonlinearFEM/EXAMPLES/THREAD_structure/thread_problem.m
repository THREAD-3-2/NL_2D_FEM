function mds = thread_problem()
%% Input definition
% mds : model data structure for automatic initialisation of the FEM model
mds = create_empty_mds();
%% geometry for a L shpae frame
% THREAD
L = 1;
% Key nodes
nodes = [1:21];
Xlist = [0 1 2 2 2 3 3 3 4 4 4 5 5 5 5.25 5.5 5.75 6 6 7 1];
Ylist = [2 0 0 1 2 0 1 2 0 1 2 0 1 2 1    2    1   0 2 1 2];
Zlist = 0*[2 0 0 1 2 0 1 2 0 1 2 0 1 2 1    2    1   0 2 1 2];
mds.geom.geom_node = [nodes', Xlist', Ylist', Zlist'];
% Lines of the geometry
lines = [1:25];
P1list = [2  1  21 3 4 4 6 7 7 7  8  9  10 9  9 10 11 12 15 16 15 17 18 18 19];
P2list = [21 21 5  4 5 7 7 8 9 10 10 10 11 11 12 13 14 15 16 17 17 18 19 20 20];
mds.geom.geom_element = [lines', P1list', P2list'];  % line 1 between point 1 and 2
% Discretisation for each line
mds.geom.discretisation = 2*ones(1,25); % line 1 discretized with 5 elements
% for user defined mesh only
mds.mesh.nodes = [];
mds.mesh.connect = [];
%% material inputs + adimentionalization
b = 20; h = 20;  % constant relative to the beam section
S = b*h;       % section area
I = b*h^3/12;  % section second moment of area
rho = 2.8e-9;  % material density
E = 70e3;     % material Young Modulus
mds.prop.poisson =0.3;  % material Poisson ratio
mds.prop.k = 1;           % shear area coeff
mds.prop.alpha = 0.1;   % Rayleigh damping coeff (C = alpha M)
% non dimensional beam material inputs + library of materials
epsi = I/S/L^2;
mds.prop.S = 1; mds.prop.I = epsi; mds.prop.rho = 1; mds.prop.E = 1/epsi;
mds.geom.geom_node(:,2:3) = mds.geom.geom_node(:,2:3)/L; % rescale geometry if adim
%% Boundary condition input (Dirichlet condition only)
mds.boundary.bc_node_list{1} = struct('node', 2,'dof', [1,2,3]); % node 1 is clamped
% mds.boundary.bc_node_list{2} = struct('node', 3,'dof', [1,2,3]); % node 1 is clamped
% mds.boundary.bc_node_list{3} = struct('node', 6,'dof', [1,2,3]); % node 1 is clamped
% mds.boundary.bc_node_list{4} = struct('node', 9,'dof', [1,2,3]); % node 1 is clamped
% mds.boundary.bc_node_list{5} = struct('node', 12,'dof', [1,2,3]); % node 1 is clamped
mds.boundary.bc_node_list{2} = struct('node', 18,'dof', [1,2,3]); % node 1 is clamped
%% Visualisation node for results
mds.visu.visu_node_list{1} = struct('node', 8 , 'dof', [1 2]);   % ux and uy of node 2 are observed
%% Forcing condition input (static and periodic)
% ponctual static forces
% mds.force.static.static_ponctual_force_node_list{1} = struct('node', [3],'dof', [2],'amplitude', [1] ); % Fx = -1 and Fy = 0 at node 2
% distrubuted static forces
mds.force.static.static_distributed_force_amplitude = [0.001];
mds.force.static.static_distributed_force_direction = [2];
% ponctual periodic forces
mds.force.periodic.periodic_ponctual_force_node_list{1} = struct('node', 8,'dof', [1],'amplitude', [0.5], 'harmonic', [1] ); % complex amplitude f = re(amp) cos + im(amp) sin
% distrubuted periodic forces
% mds.force.periodic.periodic_distributed_force_amplitude = []; % complex amplitude f = re(amp) cos + im(amp) sin
% mds.force.periodic.periodic_distributed_force_direction = [];
% mds.force.periodic.periodic_distributed_force_harmonic = [];
end