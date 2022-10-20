function mds = toy_house_problem()
%% Input definition
% mds : model data structure for automatic initialisation of the FEM model
mds = create_empty_mds();
%% geometry for a toy house frame
%   ^         (5)
%   |        / | \
%   |       /  |  \            % [l]: line numbers
%   | [4]  /   |[6]\  [5]      % (p): point Numbers
%   |     /   (6)   \
% H |    /           \
%   |  (4)----------(3)     
%   |   |     [2]    |
%   |   | [3]        | [1]
%   |   |            |
%   v  (1)          (2)
%       <------------->
%              L
%
L = 1000; % length
H = 1000; % height
% Key nodes
mds.geom.geom_node = [
    1  0   0 0;
    2  L   0 0;
    3  L    H/2 0;
    4  0    H/2 0;
    5  L/2  H 0;
    6  L/2  3*H/4 0];
% Lines of the geometry
mds.geom.geom_element = [1 2 3;
    2 3 4;
    3 4 1 ;
    4 4 5;
    5 5 3 ;
    6 5 6;];
% Discretisation for each line
mds.geom.discretisation = [3, 3, 3, 3, 3, 2];

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
mds.boundary.bc_node_list{2} = struct('node', 2,'dof', [1,2,3]); % node 2 is clamped
%% Visualisation node for results
mds.visu.visu_node_list{1} = struct('node', 5 , 'dof', [1 2]);   % ux and uy of node 2 are observed
%% Forcing condition input (static and periodic)
% ponctual static forces
mds.force.static.static_ponctual_force_node_list{1} = struct('node', [6],'dof', [2],'amplitude', [0] ); % Fx = -1 and Fy = 0 at node 2
% distrubuted static forces
% mds.force.static.static_distributed_force_amplitude = [-4];
% mds.force.static.static_distributed_force_direction = [2];
% ponctual periodic forces
% mds.force.periodic.periodic_ponctual_force_node_list{1} = struct('node', 2,'dof', [1],'amplitude', [0.1], 'harmonic', [1] ); % complex amplitude f = re(amp) cos + im(amp) sin
% distrubuted periodic forces
mds.force.periodic.periodic_distributed_force_amplitude = [0.01]; % complex amplitude f = re(amp) cos + im(amp) sin
mds.force.periodic.periodic_distributed_force_direction = [1];
mds.force.periodic.periodic_distributed_force_harmonic = [1];
end