clear
close all
clc

%% FEM model definition
% create a new NL_2D_FEM model
model = NL_2D_FEM;

L = 1;
% Key nodes
geom_node = [1 0 0 0;   % Point 1 coord (0,0)
    2 L/pi 0 0; % Point 2 coord (L/pi,0) % for force
    3  L  0 0;];% Point 3 coord (L, 0)
% Lines of the geometry
geom_element = [1 1 2    % line 1 between point 1 and 2
    2 2 3];  % line 2 between point 2 and 3
% Discretisation for each line
discretisation = [5, 10]; % line 1 discretized with 15 elements
% set geometry (optional if a mesh is directly provided)
model = model.set_geom(geom_node, geom_element, discretisation);

% set compute a mesh
model = model.auto_mesh();

%% material inputs + adimentionalization
b = 50; h = 1;  % constant relative to the beam section
S = b*h;       % section area
I = b*h^3/12;  % section second moment of area

% non dimensional parameters
poisson = 0.3;   % material Poisson ratio
k = 1;           % shear area coeff
alpha = 0.1;     % Raileight damping coeff (C = alpha M)
epsi = I/(S*L^2);% slenderness ratio
S = 1;           % nd area
I = epsi;        % nd 2nd moment
rho = 1;         % nd density
E = 1/epsi;      % nd Young mod

% set properties
model = model.set_prop(S, I, rho, E, poisson, k, alpha);

%% Boundary condition input (Dirichlet condition only)
bc_node_list{1} = struct('node', 1,'dof', [1,2,3]); % node 1 is clamped
bc_node_list{2} = struct('node', 3,'dof', [1,2,3]); % node 3 is clamped
% set boundary condition
model = model.set_boundary(bc_node_list);

%% Visualisation node Input (for results display)
visu_node_list{1} = struct('node', 2 ,...
    'dof', [1 2]);
% set visualized nodes
model = model.set_visu(visu_node_list);

% point periodic force
periodic_ponctual_force_node_list{1} = struct('node', 2,'dof', [2],'amplitude', [0.1], 'harmonic', [1] ); % complex amplitude f = re(amp) cos + im(amp) sin

% dynamic loads
model = model.set_periodic_loads('ponctual', periodic_ponctual_force_node_list);

% assemble mass matrix and force vector
model = model.initialise_matrices_and_vector();

% static solution
[qs_full, res] = model.solve_static_problem();
model.plot_static_configuration(qs_full)

% compute stresses
[strain,stress] = model.strains_and_stress_at_gauss_point(qs_full);


% modal analysis
[shape, freq] = model.linear_modal_analysis(qs_full);
% [shape, freq] = model.linear_modal_analysis();
model.plot_mode_shape([1:4],shape)

% linear analysis
H = 1;
target_mode = 1;
Omega = linspace(freq(target_mode)*0.8, freq(target_mode)*1.2,500)*2*pi;
[qp_full, bode] = model.linear_analysis(H, Omega, qs_full);
figure
subplot(2,1,1) % first harmonic amplitude; hold on
plot(Omega, bode.amp_qp_full{1}(4,:)) % u node 2
plot(Omega, bode.amp_qp_full{1}(5,:)) % v node 2
xlabel('Omega'); ylabel('Amp H1')
subplot(2,1,2) % first harmonic phase; hold on
plot(Omega, bode.phase_qp_full{1}(5,:)) % phase u node 2
plot(Omega, bode.phase_qp_full{1}(5,:)) % phase v node 2
xlabel('Omega'); ylabel('Phase H1')


%% MANLAB LAUNCHING SEQUENCE
%% MANLAB INPUTS
global U Section Diagram   % Global variables to export point from the diagram.
H = 10;          % number of harmonics for the Fourier series
type = 'autonomous'; % type of system (can be 'forced' or 'autonomous')
% type = 'forced'; % type of system (can be 'forced' or 'autonomous')
target_mode = 1; % modeto be studied in NNM (autonomous) or in FRF (forced)
angfreq = 'omega'; % 'omega' or constant value
%% Use the model to initialise MANLAB computation automatically
% MANLAB structure of parameters for equation.m
[nz, nz_aux, parameters] = model.set_MAN_parameters(H, type,  model, angfreq);
% Construct MANLAB system (matlab object)
sys = SystHBQ(nz,nz_aux,H,@equations_vector_NL_2D_FEM,@point_display,@global_display,parameters,type,'vectorial');

if strcmp(type,'autonomous')
    omega0 = (freq(target_mode)*2*pi);
    lambda0 = 0;
    idx = sys.getcoord('cos',2 ,1); % dof to be imposed amplitude
    amp = 1e-5;    % imposed amplitude
    [Z0] = model.man_initial_point(H, omega0, qs_full, amp*shape(:,target_mode));
    U0 = sys.init_U0(Z0, omega0, lambda0);
    U0 = model.solve_MAN_system_at_fixed_amplitude(U0, idx, amp, sys);
elseif strcmp(type, 'forced')
    omega0 = freq(target_mode)*2*pi*0.8;
    lambda0 = omega0; % continuation parameter initial value
    [qp_full, bode] = model.linear_analysis(H, omega0);
    [Z0] = model.man_initial_point(H, omega0, qs_full, qp_full);
    U0 = sys.init_U0(Z0, omega0, lambda0);
    U0 = model.solve_MAN_system_at_fixed_frequency(U0, omega0, sys);
end
%%% Variable displayed in the projected bifurcation diagram.
% To plot the coefficient of cos(h omega t) of variable number i with
% respect to lambda you should write as follows:
dispvars = [sys.getcoord('omega') sys.getcoord('cos',1,1);
    sys.getcoord('omega') sys.getcoord('sin',1,1);
    sys.getcoord('omega') sys.getcoord('cos',2,1);
    sys.getcoord('omega') sys.getcoord('sin',2,1)];
%% Launch of Manlab with options
Manlab('sys'       ,sys , ...
    'U0value'         ,U0, ...
    'order'           ,20, ...     % order of the series
    'ANMthreshold'    ,1e-10, ...   % threshold for the domain of validity of the series
    'Amax_max'        ,1e2, ...    % maximum value of the domain of validity of the series
    'NRthreshold'     ,1e-12, ...   % threshold for Newton-Raphson (NR) corrections
    'NRitemax'        ,50, ...     % Maximum number of iteration of NR algorithm
    'NRstart'         ,0, ...      % NR corrections for the user-defined starting point [on]/off
    'NRmethod'        ,0, ...      % NR corrections on/[off]
    'BifDetection'    ,1, ...      % Detection of bifurcation [on]/off
    'PointDisplay'    ,0, ...      % Point display [on]/off
    'GlobalDisplay'   ,0, ...      % Global display [on]/off
    'StabilityCheck'  ,0, ...      % Stability computation on/[off]
    'StabTol'         ,1e-6, ...   % Stability tolerance
    'displayvariables',dispvars);     % MANLAB run