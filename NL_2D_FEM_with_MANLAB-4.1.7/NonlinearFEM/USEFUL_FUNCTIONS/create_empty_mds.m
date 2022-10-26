function mds = create_empty_mds()
% mds : model data structure for automatic initialisation of the FEM model
mds.geom.geom_node = [ ];
mds.geom.geom_element = [];
mds.geom.discretisation = [];
%% for user defined mesh only
mds.mesh.nodes = [];
mds.mesh.connect = [];
%% material inputs + adimentionalization
mds.prop.S = []; 
mds.prop.I = []; 
mds.prop.k = [];           % shear area coeff
mds.prop.rho = [];
mds.prop.E = [];
mds.prop.poisson =[];  % material Poisson ratio
mds.prop.alpha = [];   % Raileight damping coeff (C = alpha M)
%% Boundary condition input (Dirichlet condition only)
mds.boundary.bc_node_list = {}; % node 1 is clamped
%% Visualisation node for results
mds.visu.visu_node_list = {};   % ux and uy of node 2 are observed
%% Forcing condition input (static and periodic)
% ponctual static forces
mds.force.static.static_ponctual_force_node_list = {};
% distrubuted static forces
mds.force.static.static_distributed_force_amplitude = [];
mds.force.static.static_distributed_force_direction = [];
% ponctual periodic forces
mds.force.periodic.periodic_ponctual_force_node_list = {};
% distrubuted periodic forces
mds.force.periodic.periodic_distributed_force_amplitude = []; % complex amplitude f = re(amp) cos + im(amp) sin
mds.force.periodic.periodic_distributed_force_direction = [];
mds.force.periodic.periodic_distributed_force_harmonic = [];
