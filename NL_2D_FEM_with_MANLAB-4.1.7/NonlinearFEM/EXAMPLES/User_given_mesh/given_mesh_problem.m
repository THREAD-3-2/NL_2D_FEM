function mds = given_mesh_problem(N)
%% Input definition
% mds : model data structure for automatic initialisation of the FEM model
mds = create_empty_mds();
%% geometry for a horizontal beam
%       [N-1 elements]             % line numbers
%    (1)--...(n)...--(N)     % point Numbers
%     <-------------->
%            L
%  
L=1000;
mds.geom.geom_node = [];
mds.geom.geom_element = [];
mds.geom.discretisation = [];
% user defined mesh
mds.mesh.nodes = [[1:N]',L*[0:N-1]'/(N-1),zeros(N,1)]; % nodes coordiates
mds.mesh.connect = [[1:N-1]',[1:N-1]',[2:N]' ]; % element list
%% material inputs + adimentionalization
b = 10; h = 10;  % constant relative to the beam section
S = b*h;       % section area
I = b*h^3/12;  % section second moment of area
rho = 7.8e-9;  % material density
E = 210e3;     % material Young Modulus
mds.prop.poisson =0.3;  % material Poisson ratio
mds.prop.k = 1;           % shear area coeff
mds.prop.alpha = 0.1;   % Raileight damping coeff (C = alpha M)
% non dimensional beam material inputs + library of materials
epsi = I/(S*L^2);
mds.prop.S = 1; mds.prop.I = epsi; mds.prop.rho = 1; mds.prop.E = 1/epsi;
mds.mesh.nodes(:,2:3) = mds.mesh.nodes(:,2:3)/L ;%resacle mesh
%% Boundary condition input (Dirichlet condition only)
mds.boundary.bc_node_list{1} = struct('node', 1,'dof', [1,2,3]); % node 1 is clamped
%% Visualisation node for results
mds.visu.visu_node_list{1} = struct('node', N , 'dof', [1 2]);   % ux and uy of node 2 are observed
%% Forcing condition input (static and periodic)
% ponctual static forces
mds.force.static.static_ponctual_force_node_list{1} = struct('node', N,'dof', [1 2],'amplitude', [0 0] ); % Fx = -1 and Fy = 0 at node 2
% distrubuted static forces
mds.force.static.static_distributed_force_amplitude = [10]; % distributed torque
mds.force.static.static_distributed_force_direction = [3]; % rotation direction (=torque)
% ponctual periodic forces
mds.force.periodic.periodic_ponctual_force_node_list{1} = struct('node', N,'dof', [2],'amplitude', [0], 'harmonic', [1] ); % complex amplitude f = re(amp) cos + im(amp) sin
% distrubuted periodic forces
mds.force.periodic.periodic_distributed_force_amplitude = [0 0.1+0.1*1i]; % complex amplitude f = re(amp) cos + im(amp) sin
mds.force.periodic.periodic_distributed_force_direction = [1  2];
mds.force.periodic.periodic_distributed_force_harmonic = [1  1];
end