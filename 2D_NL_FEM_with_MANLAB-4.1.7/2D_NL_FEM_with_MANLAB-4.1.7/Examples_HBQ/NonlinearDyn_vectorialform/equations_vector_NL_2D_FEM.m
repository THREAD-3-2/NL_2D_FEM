function [Rf,dRa,Forcing] = equations_vector_NL_2D_FEM(sys, t, zf, dzf, d2zf)
% Equations of the system of the form 
% Rf(zf) = C + L0(zf) + lambda L1(zf) + D0(dzf) + lambda D1(dzf) + DD(d2zf) + Q(zf,zf) + f(zf).

sizing = size(zf,2);

if strcmp(sys.parameters.type,'autonomous')
    NNM=1;
    FRF=0;
elseif strcmp(sys.parameters.type,'forced')
    NNM=0;
    FRF=1;
else
    
end
model = sys.parameters.model;
%keyboard
R     = zeros(sys.nz,sizing);      % Main residue
Ra    = zeros(sys.nz_aux,sizing);  % Auxiliary residue
dRa   = zeros(sys.nz_aux,sizing);  % Differential form of non-quadratic part of the auxiliary residue
active_dof = sys.parameters.model.boundary.active_dof;

%%% parameters of the system
number_nodes = sys.parameters.model.mesh.number_nodes;
number_elements = sys.parameters.model.mesh.number_elements;
n_var = sys.nz/2; % location of first variable u

%% Variables : 
u      = zf(1:n_var,:);           du     = dzf(1:n_var,:);             % Main variable: u = position
v      = zf(n_var + 1:2*n_var,:); dv     = dzf(n_var + 1:2*n_var,:);   % Main variable: v = velocity
aux    = zf(2*n_var + 1:end-1,:); daux   = dzf(2*n_var + 1:end,:);     % Auxiliary variables
lambda = zf(end,:);                                                  % Continuation parameter
% auxiliary variable used for diffeential equation
theta  = aux( 3*number_elements + 1: 4*number_elements,:);  dtheta = daux(3*number_elements + 1:4*number_elements,:);
c      = aux( 4*number_elements + 1: 5*number_elements,:);      dc = daux(4*number_elements + 1:5*number_elements,:);
s      = aux( 5*number_elements + 1: 6*number_elements,:);      ds = daux(5*number_elements + 1:6*number_elements,:);

% constructing displacement vector q from Main variable u:
q_full = zeros(3*number_nodes,sizing);
q_full(active_dof,:) = u;
% constructing velicity vector v from Main variable v:
v_full = zeros(3*number_nodes,sizing); % needed to calculate the residual
v_full(active_dof,:) = v;
% constructing acceleration vector dv from Main variable dv:
dv_full = zeros(3*number_nodes,sizing); % needed to calculate the residual
dv_full(active_dof,:) = dv;

% matrices
M_full =  model.matrices.mass;
C_full =  model.matrices.damping;
%static forces
f_stat_full = model.vectors.static_forces;
% internal force vector (constructed from aux var)
f_int_full = model.assemble_constant_vector('internal_force_vector_from_aux_var', aux);

% equation of motion
if NNM % NNM equation
    Residual = f_stat_full -f_int_full - M_full*dv_full - lambda.*v_full;
else % forced response equation
    Residual = f_stat_full -f_int_full - M_full*dv_full - C_full*v_full;
end
% auxiliary residual
aux_expression = model.auxiliary_variable_expression(q_full, aux);
Ra = aux-aux_expression;
% Differential equation for cosine and sine
dRa(4*number_elements + 1: 5*number_elements,:) = dc + s.*dtheta;
dRa(5*number_elements + 1: 6*number_elements,:) = ds - c.*dtheta;
%% Final Residal
% Equations of the main system
R(1:n_var,:) = v - du;
R(n_var + 1:2*n_var,:) = Residual(active_dof,:);
% Concatenation of the two residues
Rf = [ R ; Ra ];

%% Forcing terms
Forcing = zeros(2*sys.H+1,sys.nz_tot); % DO NOT CHANGE this line.
% Should be written as if the forcing angular frequency value is 1
% i.e. the forcing period is 2*pi
if FRF
    f_dyn_full = model.vectors.periodic_forces;
    f_dyn = f_dyn_full(active_dof,:);
    for h=1:size(f_dyn_full,2)
        Forcing(:,n_var + 1:2*n_var) = Forcing(:,n_var + 1:2*n_var) + cos(h*t')*real(f_dyn(:,h)') + sin(h*t')* imag(f_dyn(:,h)');
    end
end
end