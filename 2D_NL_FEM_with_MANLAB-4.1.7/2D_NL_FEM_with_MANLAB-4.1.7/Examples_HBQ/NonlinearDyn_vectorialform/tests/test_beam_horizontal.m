clear
close all
clc

%% Input definition
%% geometry for a vertical beaam
L = 1000; % beam length
% key nodes of the geometry
geom_node = [1 0 0 0;  % Point 1 coord (0,0)
             2 L 0 0;];% Point 2 coord (0,1000)         
% line of the geometry
geom_element = [1 1 2];% line 1 between point 1 and 2
discretisation = [5]; % line 1 discretized with 5 elements                  
%% material inputs
b = 10; h = 10;
S = b*h; I = b*h^3/12; rho = 7.8e-9; E = 210e3; poisson =0.3; k=5/6;
alpha = 0.1; % Raileight damping coeff C = alpha M
% non dimensional beam material
epsi = I/(S*L^2);
S = 1; I = epsi; rho = 1; E = 1/epsi; poisson =0.3; k=1;
geom_node(:,2:3) = geom_node(:,2:3)/L; % adimentionalize geometry
%% Boundary condition input (Dirichlet condition only)
bc_node_list{1} = struct('node', 1,'dof', [1,2,3]);
% bc_node_list{2} = struct('node_number',2,'prescribed_dof',[1,2,3]); % uncomment for clamped clamped beam
%% Visualisation node Input (for results display)
visu_node_list{1} = struct('node', 2 ,...
                            'dof', [1 2]);
%% Forcing Definition condition input (Dirichlet condition only)
% ponctual static forces
static_ponctual_force_node_list{1} = struct('node',2,...
                                     'dof', [1 2],...
                                     'amplitude', [-1 0] );
% distrubuted static forces
static_distributed_force_amplitude = [0  0];
static_distributed_force_direction = [1  2];
% ponctual periodic forces
periodic_ponctual_force_node_list{1} = struct('node',2,...
                                        'dof', [2],...
                                        'amplitude', [0.001],... % complex amplitude re cos + im sin 
                                        'harmonic', [1] );
% distrubuted periodic forces
periodic_distributed_force_amplitude = [0 0]; % complex amplitude
periodic_distributed_force_direction = [1  2];
periodic_distributed_force_harmonic = [1  1];
%% FEM model construction
model = NL_2D_FEM; % create a new NL_2D_FEM model
% set geometry (optional if a mesh is directly provided)
model = model.set_geom(geom_node, geom_element, discretisation);
% set visualized nodes
model = model.set_visu(visu_node_list); 
% set mesh 
model = model.auto_mesh(); % auto mesh if geometry is given
% model = model.set_mesh(nodes, connect; % user defined mesh otherwise: needs nodes coordiantes and connectivity table
% set properties
model = model.set_prop(S, I, rho, E, poisson, k, alpha); % set propoerties (must be the same for all elements)
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

[qs_full, res] = model.solve_static_problem();
q0_full = model.vectors.null_vector;
fig = figure(99);
model.plot_deformed_mesh(q0_full, fig, '-k')
model.plot_deformed_mesh(qs_full, fig, '-b')
[strain,stress] = model.strains_and_stress_at_gauss_point(qs_full);
keyboard
% modal analysis
[shape, freq] = model.linear_modal_analysis(qs_full);
% analytical frequencies
fax = sqrt(E/rho)*pi/2/L*[1, 3, 5]/2/pi
fflex = sqrt(E*I/rho/S/L^4)*[1.87^2, 4.69^2 7.85^2 10.99^2]/2/pi

% [shape, freq] = model.linear_modal_analysis();
fig3 = figure(98);
model.plot_deformed_mesh(q0_full, fig3, '--k')
model.plot_deformed_mesh(qs_full, fig3, '--b')
model.plot_deformed_mesh(shape(:,1), fig3, '-b')
model.plot_deformed_mesh(shape(:,2), fig3, '-g')
model.plot_deformed_mesh(shape(:,3), fig3, '-g')
% buckling analysis
% [frq_list, lambda, N_list, q_list] = model.buckling_analysis();
% figure; plot(lambda, frq_list(1,:))
% figure; hold on; plot(lambda, q_list(4,:)); plot(lambda, q_list(5,:))
% fig44 = figure(44)
% hold on
% model.plot_deformed_mesh(q_list(:,1), fig44, '-b')
% model.plot_deformed_mesh(q_list(:,floor(end/2)), fig44, '-k')
% model.plot_deformed_mesh(q_list(:,end), fig44, '-g')
% 
% K_tan = model.assemble_constant_matrix('stiffness_at_qs', qs_full);
% K_lin = model.assemble_constant_matrix('stiffness_at_qs', 0*qs_full);
% 
% keyboard
%%

global U Section Diagram   % Global variables to export point from the diagram.

% Path of the SRC file.
addpath(genpath('..\..\NonlinearDyn_vectorialform'));
addpath(genpath('..\..\..\SRC'));

%% Parameters of the system

neq = length(model.boundary.active_dof); %% number of variables: [u,w,theta] per node
neq_aux = 12*model.mesh.number_elements; %% number of auxiliary variables: in this case, 12 per element

%% initialization of the system
H = 15;             % number of harmonics used to compute the solution-branch
nz = 2*neq;         % number of main equations of the system of the differential-algebraic system (DAE)
nz_aux = neq_aux;   % number of auxiliary equations of the system of the DAE


Omega = linspace(freq(1)*0.8, freq(1)*1.2, 500)*2*pi;
% linear analysis
[qp_full, bode] = model.linear_analysis(H, Omega, qs_full);
figure
subplot(2,1,1)
hold on
plot(Omega, bode.amp_qp_full{1}(4,:)) % u
plot(Omega, bode.amp_qp_full{1}(5,:)) % v
plot(Omega, bode.amp_qp_full{1}(6,:)) % theta
xlabel('Omega'); ylabel('Amp H1')
subplot(2,1,2)
hold on
plot(Omega, bode.phase_qp_full{1}(4,:)) % u
plot(Omega, bode.phase_qp_full{1}(5,:)) % v
plot(Omega, bode.phase_qp_full{1}(6,:)) % theta
xlabel('Omega'); ylabel('Phase H1')
% figure
% subplot(2,1,1)
% hold on
% plot(Omega, bode.amp_qp_full{2}(4,:)) % u
% plot(Omega, bode.amp_qp_full{2}(5,:)) % v
% plot(Omega, bode.amp_qp_full{2}(6,:)) % theta
% xlabel('Omega'); ylabel('Amp H2')
% subplot(2,1,2)
% hold on
% plot(Omega, bode.phase_qp_full{2}(4,:)) % u
% plot(Omega, bode.phase_qp_full{2}(5,:)) % v
% plot(Omega, bode.phase_qp_full{2}(6,:)) % theta
% xlabel('Omega'); ylabel('Phase H2')

% MANLAB structure of parameters for equation.m
parameters.H = H;
parameters.angfreq = 'omega'; % if the system is forced at a fixed angular frequency, parameters.angfreq
                              % contains its value. Otherwise, parameters.angfreq = 'omega'.
type = 'autonomous'; % type of system (can be 'forced' or 'autonomous')
%type = 'forced'; % type of system (can be 'forced' or 'autonomous')
writing = 'vectorial'; % way the equations have been written (can be 'standard' or 'vectorial')
parameters.type = type;
parameters.model = model;
sys = SystHBQ(nz,nz_aux,H,@equations_vector_NL_2D_FEM,@point_display,@global_display,parameters,type,writing);

omega0 = (freq(1)*2*pi);
% lambda0 = omega0; % continuation parameter initial value
lambda0 = 0; % continuation parameter initial value
epsilon = 0.01; % amplitude of the linear mode

%% starting point
[qp0_full, bode] = model.linear_analysis(H, omega0);
qpH1c = bode.qcos(:,1);
qpH1s = bode.qsin(:,1);

% [Z0] = model.man_initial_point(H, omega0, qs_full, -qp0_full);
[Z0] = model.man_initial_point(H, omega0, qs_full, epsilon*shape(:,1));

fig55=figure(55);
hold on
model.plot_deformed_mesh(qs_full, fig55, '--k')
model.plot_deformed_mesh([0;0;0;Z0(2,1:length(model.boundary.active_dof))'], fig55, '--ro')
model.plot_deformed_mesh([0;0;0;Z0(2+H,1:length(model.boundary.active_dof))'], fig55, '--bo')
U0z = sys.init_U0(Z0, omega0, lambda0);

% % solve system at stating frequency 
% options = optimoptions('fsolve','SpecifyObjectiveGradient',true,'display','iter','maxiter',30,'FunctionTolerance',1e-8,'StepTolerance',1e-12, 'Algorith','levenberg-marquardt');% trust-region-dogleg
% U0 = fsolve(@(U)man_residual_fixed_frequency(U, omega0, sys), U0z, options);
% U01 = fsolve(@(U)man_residual_fixed_frequency(U, omega0, sys), U0, options);
% Res = man_residual_fixed_frequency(U0, omega0, sys);
% fprintf(1,'HBM system at initial frequency solved, residual %2.2e \n', norm(Res))
% [Ztot,omega,lambda,omega2,lambdaomega] = get_Ztot(sys,U01);
% model.plot_deformed_mesh([0;0;0;Ztot(2,1:length(model.boundary.active_dof))'], fig55, '-m*')
% model.plot_deformed_mesh([0;0;0;Ztot(2+H,1:length(model.boundary.active_dof))'], fig55, '-c*')

amp = 1e-2;
idx = sys.getcoord('cos',2,1);
% solve system at stating amplitude 
options = optimoptions('fsolve','SpecifyObjectiveGradient',true,'display','iter','maxiter',30,'FunctionTolerance',1e-8,'StepTolerance',1e-12, 'Algorith','levenberg-marquardt');% trust-region-dogleg
U0 = fsolve(@(U)man_residual_fixed_amplitude(U, idx, amp, sys), U0z, options);
U01 = fsolve(@(U)man_residual_fixed_amplitude(U, idx, amp, sys), U0, options);
Res = man_residual_fixed_amplitude(U0, idx, amp, sys);
fprintf(1,'HBM system at initial amplitude solved, residual %2.2e \n', norm(Res))
[Ztot,omega,lambda,omega2,lambdaomega] = get_Ztot(sys,U01);
model.plot_deformed_mesh([0;0;0;Ztot(2,1:length(model.boundary.active_dof))'], fig55, '-m*')
model.plot_deformed_mesh([0;0;0;Ztot(2+H,1:length(model.boundary.active_dof))'], fig55, '-c*')


% model.plot_deformed_mesh(qpH1c,fig55,'-*')
% model.plot_deformed_mesh(qpH1s,fig55,'-*')

U0 = sys.init_U0(Ztot, omega0, lambda0);
max(abs(U0-U01))

% model.plot_deformed_mesh([0;0;0;Ztot(2,1:length(model.boundary.active_dof))'], fig, '-r')
% model.plot_deformed_mesh([0;0;0;Ztot(2+H,1:length(model.boundary.active_dof))'], fig, '-r')
% model.plot_deformed_mesh(qpH1c, fig, '-b')
% model.plot_deformed_mesh(qpH1s, fig, '-k')

% u_full = real(qp0_full(:,1))-imag(qp0_full(:,1));
% v_full = -omega0*u_full;
% u = u_full(model.boundary.active_dof);
% v = v_full(model.boundary.active_dof);
% aux = model.man_auxiliary_variables_vector(u_full);
% zf=[u; v; aux; omega0];
% Res = sys.equations(sys, 0, zf, 0*zf, 0*zf)

%%% Phase initialization for nonlinear modes 
% sys.zi_phase = nz/2 - 1;
% sys.zi_phase = n_var2;
%%% Variable displayed in the projected bifurcation diagram.
% To plot the coefficient of cos(h omega t) of variable number i with
% respect to lambda you should write as follows:
dispvars = [sys.getcoord('omega') sys.getcoord('cos',1,1);
            sys.getcoord('omega') sys.getcoord('sin',1,1);
            sys.getcoord('omega') sys.getcoord('cos',2,1);
            sys.getcoord('omega') sys.getcoord('sin',2,1)];
%             sys.getcoord('omega') sys.getcoord('cos',dof_info.obs_dof,3);
%             sys.getcoord('omega') sys.getcoord('sin',dof_info.obs_dof,3)];
%             sys.getcoord('omega') sys.getcoord('cos',dof_info.obs_dof,5);
%             sys.getcoord('omega') sys.getcoord('sin',dof_info.obs_dof,5)];

%% Launch of Manlab with options
Manlab('sys'       ,sys , ...
    'U0value'         ,U0, ...
    'order'           ,20, ...     % order of the series
     'ANMthreshold'    ,1e-10, ...   % threshold for the domain of validity of the series
     'Amax_max'        ,1e2, ...    % maximum value of the domain of validity of the series
     'NRthreshold'     ,1e-12, ...   % threshold for Newton-Raphson (NR) corrections
     'NRitemax'        ,50, ...     % Maximum number of iteration of NR algorithm
     'NRstart'         ,1, ...      % NR corrections for the user-defined starting point [on]/off
     'NRmethod'        ,0, ...      % NR corrections on/[off]
     'BifDetection'    ,1, ...      % Detection of bifurcation [on]/off
     'PointDisplay'    ,0, ...      % Point display [on]/off
     'GlobalDisplay'   ,0, ...      % Global display [on]/off
     'StabilityCheck'  ,0, ...      % Stability computation on/[off]
     'StabTol'         ,1e-6, ...   % Stability tolerance
    'displayvariables',dispvars);     % MANLAB run


%% Optional arguments of Manlab launching function with their default values :
%          'order'           ,20, ...     % order of the series
%          'ANMthreshold'    ,1e-6, ...   % threshold for the domain of validity of the series
%          'Amax_max'        ,1e6, ...    % maximum value of the domain of validity of the series
%          'NRthreshold'     ,2e-5, ...   % threshold for Newton-Raphson (NR) corrections
%          'NRitemax'        ,10, ...     % Maximum number of iteration of NR algorithm
%          'NRstart'         ,1, ...      % NR corrections for the user-defined starting point [on]/off
%          'NRmethod'        ,0, ...      % NR corrections on/[off]
%          'BifDetection'    ,1, ...      % Detection of bifurcation [on]/off
%          'PointDisplay'    ,1, ...      % Point display [on]/off
%          'GlobalDisplay'   ,1, ...      % Global display [on]/off
%          'StabilityCheck'  ,0, ...      % Stability computation on/[off]
%          'StabTol'         ,1e-6, ...   % Stability tolerance
%
%               The stability requires to write the system of equations in
%               the explicit ODE form X' = f(X).
%               The residue function is then R(X) = f(X) = 0.
%               In all other cases, it gives wrong results.
