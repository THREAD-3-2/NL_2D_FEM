function mds = L_shape_problem()
%% Input definition
% mds : model data structure for automatic initialisation of the FEM model
mds = create_empty_mds();
%% geometry for a L shpae frame
%             [2]             % line numbers
%   ^ (5)-----(6)------------(7)-----(8)     % point Numbers
%   |          |              |
%   |          |              |
% H |          | [1]          |
%   |          |              |
%   |          |              |
%   v (1)-----(2)------------(3)-----(4)
%       <----------------------------->
%                        L
%  
L = 1000; % length
H = 2000; % height
% Key nodes
mds.geom.geom_node = [1 0 0 0;  % Point 1 coord (0,0)
                      2 L/4 0 0;  % Point 2 coord (0,H)
                      3 3*L/4 0 0;
                      4 L 0 0
                      5 0 H 0;  % Point 1 coord (0,0)
                      6 L/4 H 0;  % Point 2 coord (0,H)
                      7 3*L/4 H 0;
                      8 L H 0;];% Point 3 coord (L,H)
% Lines of the geometry
mds.geom.geom_element = [1 1 2
                         2 2 3
                         3 3 4
                         4 5 6
                         5 6 7
                         6 7 8
                         7 2 6 
                         8 3 7];  % line 1 between point 1 and 2
% Discretisation for each line
mds.geom.discretisation = [2 4 2 2 4 2 6 6]; % line 1 discretized with 5 elements
% for user defined mesh only
mds.mesh.nodes = [];
mds.mesh.connect = [];
%% material inputs + adimentionalization
b = 10; h = 10;  % constant relative to the beam section
S = b*h;       % section area
I = b*h^3/12;  % section second moment of area
rho = 7.8e-9;  % material density
E = 210e3;     % material Young Modulus
mds.prop.poisson =0.3;  % material Poisson ratio
mds.prop.k = 1;           % shear area coeff
mds.prop.alpha = 0.1;   % Rayleigh damping coeff (C = alpha M)
% non dimensional beam material inputs + library of materials
epsi = I/(S*L^2);
mds.prop.S = 1; mds.prop.I = epsi; mds.prop.rho = 1; mds.prop.E = 1/epsi;
mds.geom.geom_node(:,2:3) = mds.geom.geom_node(:,2:3)/L; % rescale geometry if adim
%% Boundary condition input (Dirichlet condition only)
mds.boundary.bc_node_list{1} = struct('node', 1,'dof', [1,2,3]); % node 1 is clamped
mds.boundary.bc_node_list{2} = struct('node', 4,'dof', [1,2,3]); % node 1 is clamped
mds.boundary.bc_node_list{3} = struct('node', 5,'dof', [1,2,3]); % node 1 is clamped
mds.boundary.bc_node_list{4} = struct('node', 8,'dof', [1,2,3]); % node 1 is clamped
%% Visualisation node for results
mds.visu.visu_node_list{1} = struct('node', 2 , 'dof', [1 2]);   % ux and uy of node 2 are observed
%% Forcing condition input (static and periodic)
% ponctual static forces
mds.force.static.static_ponctual_force_node_list{1} = struct('node', [],'dof', [],'amplitude', [] ); % Fx = -1 and Fy = 0 at node 2
% distrubuted static forces
mds.force.static.static_distributed_force_amplitude = [1];
mds.force.static.static_distributed_force_direction = [1];
% ponctual periodic forces
% mds.force.periodic.periodic_ponctual_force_node_list{1} = struct('node', [],'dof', [1],'amplitude', [0.5], 'harmonic', [1] ); % complex amplitude f = re(amp) cos + im(amp) sin
% distrubuted periodic forces
% mds.force.periodic.periodic_distributed_force_amplitude = []; % complex amplitude f = re(amp) cos + im(amp) sin
% mds.force.periodic.periodic_distributed_force_direction = [];
% mds.force.periodic.periodic_distributed_force_harmonic = [];
end