clear
close all
clc
% dont forget to set the path !!!
%% FEM model definition
% load the model data structure
mds = H_shape_problem();
% create a new NL_2D_FEM model
model = NL_2D_FEM(); 
% set up the model using the data structure
model = model.set_model_from_mds(mds);
%% MANLAB LAUNCHING SEQUENCE
global U Section Diagram   % Global variables to export point from the diagram in GUI
%% MANLAB INPUTS
H = 10;          % number of harmonics for the Fourier series   
% type = 'autonomous'; % type of system (can be 'forced' or 'autonomous')
type = 'forced'; % type of system (can be 'forced' or 'autonomous')
target_mode = 1; % modeto be studied in NNM (autonomous) or in FRF (forced)
angfreq = 'omega'; % 'omega' or constant value 
%% Use the model to initialise MANLAB computation automatically
% MANLAB structure of parameters for equation.m
[nz, nz_aux, parameters] = model.set_MAN_parameters(H, type,  model, angfreq);
% Construct MANLAB system (matlab object)
sys = SystHBQ(nz,nz_aux,H,@equations_vector_NL_2D_FEM,@point_display,@global_display,parameters,type,'vectorial');
% compute static equilibrium, modal analysis and the MANLAB starting point
[U0, omega0, lambda0] = model.initialise_MAN_computation(sys, type, target_mode);
%%% Variable displayed in the projected bifurcation diagram (GUI).
dispvars = [sys.getcoord('omega') sys.getcoord('cos',1,1);  % -> sys.getcoord('sin', dof_index, harmonic_number)
            sys.getcoord('omega') sys.getcoord('sin',1,1);
            sys.getcoord('omega') sys.getcoord('cos',2,1);
            sys.getcoord('omega') sys.getcoord('sin',2,1)];
%% Call the MANLAB GUI for numerical continuation of the problem w.r.t. Omega
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
