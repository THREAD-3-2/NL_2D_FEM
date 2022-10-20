clear
close all
clc

% dont forget to set the path !!!
%% FEM model definition
% load the model data structure
mds = H_shape_problem();
% create a new NL_2D_FEM model
model = NL_2D_FEM; 
% set up the model using the data structure
model = model.set_model_from_mds(mds);
%% Manual Use of the class NL2D_FEM
% assemble mass matrix and force vector
model = model.initialise_matrices_and_vector();
% static solution
fig_static = figure; fig_static.Name='Static displacement';
[qs_full, res] = model.solve_static_problem();
model.plot_deformed_mesh(model.vectors.null_vector, fig_static, '--k') % undeformed mesh
model.plot_deformed_mesh(qs_full, fig_static, '-b') % static deformation
[static_strain, static_stress] = model.strains_and_stress_at_gauss_point(qs_full); % static stain

model.plot_deformed_mesh_with_stress(qs_full, static_stress, fig_static, '-b') % static deformation

% modal analysis around static configuration
[shape, freq] = model.linear_modal_analysis(qs_full); % use tangent stiffnes matrices
fig_mode = figure(98); fig_mode.Name='Mode shape';
model.plot_deformed_mesh(model.vectors.null_vector, fig_mode, '--k')  % undeformed mesh
model.plot_deformed_mesh(shape(:,1), fig_mode, '-b') % mode shape display
model.plot_deformed_mesh(shape(:,2), fig_mode, '-g')
model.plot_deformed_mesh(shape(:,3), fig_mode, '-r')

%% MANLAB LAUNCHING SEQUENCE

global U Section Diagram   % Global variables to export point from the diagram.
%% MANLAB INPUTS
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

%% starting point
if strcmp(type,'forced')
    omega0 = freq(target_mode)*2*pi*0.8;
    lambda0 = omega0; % continuation parameter initial value
    [qp_full, bode] = model.linear_analysis(H, omega0);    
    [Z0] = model.man_initial_point(H, omega0, qs_full, qp_full);
    U0 = sys.init_U0(Z0, omega0, lambda0);
    U0 = model.solve_MAN_system_at_fixed_frequency(U0, omega0, sys);
elseif strcmp(type,'autonomous')
    omega0 = (freq(target_mode)*2*pi);
    lambda0 = 0;
    idx = sys.getcoord('cos',2,1); % dof to be imposed amplitude
    amp = 1e-5;    % imposed amplitude
    [Z0] = model.man_initial_point(H, omega0, qs_full, amp*shape(:,target_mode));
    U0 = sys.init_U0(Z0, omega0, lambda0);
    U0 = model.solve_MAN_system_at_fixed_amplitude(U0, idx, amp, sys);    
end
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
