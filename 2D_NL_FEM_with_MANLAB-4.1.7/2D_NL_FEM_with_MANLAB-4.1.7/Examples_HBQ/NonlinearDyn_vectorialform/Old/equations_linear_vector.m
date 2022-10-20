function [Rf,dRa,Forcing] = equations_linear_vector(sys, t, zf, dzf, d2zf)
% Equations of the system of the form 
% Rf(zf) = C + L0(zf) + lambda L1(zf) + D0(dzf) + lambda D1(dzf) + DD(d2zf) + Q(zf,zf) + f(zf).

sizing = size(zf,2);

%keyboard
R     = zeros(sys.nz,sizing);      % Main residue
Ra    = zeros(sys.nz_aux,sizing);  % Auxiliary residue
dRa   = zeros(sys.nz_aux,sizing);  % Differential form of non-quadratic part of the auxiliary residue

%%% parameters of the system
number_nodes = sys.parameters.mesh.number_nodes;
number_elements = sys.parameters.mesh.number_elements;
n_var = sys.nz/2; % location of first variable u

%index = sys.parameters.element_geometry.element_index;
active_dof = sys.parameters.dof_info.active;
prescribed_dof = sys.parameters.dof_info.bc.prescribed_dof;

E = sys.parameters.beam_params.modulus;
A = sys.parameters.beam_params.area;
I = sys.parameters.beam_params.I;
G = sys.parameters.beam_params.shear;
k = sys.parameters.beam_params.k; % shear coefficient
% 
% L_element = sys.parameters.element_geometry.L_element;
% theta_element = sys.parameters.element_geometry.theta_element;

%% Variables : 
u      = zf(1:n_var,:);           du     = dzf(1:n_var,:);             % Main variable: u = position
v      = zf(n_var + 1:2*n_var,:); dv     = dzf(n_var + 1:2*n_var,:);   % Main variable: v = velocity
aux    = zf(2*n_var + 1:end-1,:); daux   = dzf(2*n_var + 1:end,:);     % Auxiliary variables
lambda = zf(end,:);                                                  % Continuation parameter

% decomposition of aux:
up     = aux(              1:   number_elements,:);
wp     = aux(   number_elements + 1: 2*number_elements,:);
thetap = aux( 2*number_elements + 1: 3*number_elements,:);
theta  = aux( 3*number_elements + 1: 4*number_elements,:);  
c      = aux( 4*number_elements + 1: 5*number_elements,:);      
s      = aux( 5*number_elements + 1: 6*number_elements,:);      
eps    = aux( 6*number_elements + 1: 7*number_elements,:);
gam    = aux( 7*number_elements + 1: 8*number_elements,:);
Fx     = aux( 8*number_elements + 1: 9*number_elements,:);
Fy     = aux( 9*number_elements + 1:10*number_elements,:);
M      = aux(10*number_elements + 1:11*number_elements,:);
T2     = aux(11*number_elements + 1:12*number_elements,:);

% constructing q from Main variable u:
q_global = zeros(3*number_nodes,sizing);
q_global(prescribed_dof,:) = 0; %%%% generalize (see Cochelin line 59)
q_global(active_dof,:) = u;

f_int_global = zeros(3*number_nodes,sizing);
M_global =  sys.parameters.M_glob;

for i = 1:number_elements
    nodeA = sys.parameters.mesh.connect(i,2);
    nodeB = sys.parameters.mesh.connect(i,3);
    
    index_global_A = 3*(nodeA - 1) + (1:3);
    index_global_B = 3*(nodeB - 1) + (1:3);
    
    index = [index_global_A index_global_B]';
    
    [L_element(i), theta_element(i)] = elementOrientation(sys.parameters.mesh,i);
    rot_matrix = transformationMatrix(theta_element(i));
    q_element = rot_matrix'*q_global(index,:); % rotating q_element into local frame
    
    f_int_elem = [-Fx(i,:); -Fy(i,:); -M(i,:) + T2(i,:)*L_element(i)/2; Fx(i,:); Fy(i,:); M(i,:) + T2(i,:)*L_element(i)/2];
    f_int_global(index,:) = f_int_global(index,:) + rot_matrix*f_int_elem;

    Ra(i,:)                 = up(i,:)     - ((q_element(4,:) - q_element(1,:))/L_element(i));   
    Ra(i +   number_elements,:) = wp(i,:)     - ((q_element(5,:) - q_element(2,:))/L_element(i));
    Ra(i + 2*number_elements,:) = thetap(i,:) - ((q_element(6,:) - q_element(3,:))/L_element(i));
    Ra(i + 3*number_elements,:) = theta(i,:)  - ((q_element(6,:) + q_element(3,:))/2);    
    
end
N = E*A*eps;
T = k*G*A*gam;

% Remaining auxiliary variable residuals
Ra( 4*number_elements + 1: 5*number_elements,:)  = c - 1;                 
Ra( 5*number_elements + 1: 6*number_elements,:)  = s - theta;                 
Ra( 6*number_elements + 1: 7*number_elements,:)  = eps -  up;
Ra( 7*number_elements + 1: 8*number_elements,:)  = gam - wp + s;
Ra( 8*number_elements + 1: 9*number_elements,:)  = Fx - N ;
Ra( 9*number_elements + 1: 10*number_elements,:) = Fy - T;
Ra(10*number_elements + 1: 11*number_elements,:) = M - E*I*thetap;
Ra(11*number_elements + 1: 12*number_elements,:) = T2 - N.*gam + (1+eps).*T;

C_global = sys.parameters.C_glob;

v_global = zeros(3*number_nodes,sizing); % needed to calculate the residual
v_global(active_dof,:) = v;

dv_global = zeros(3*number_nodes,sizing); % needed to calculate the residual
dv_global(active_dof,:) = dv;

Residual = -f_int_global - M_global*dv_global - C_global*v_global;
%% Residues
% Equations of the main system
R(1:n_var,:) = v - du;
R(n_var + 1:2*n_var,:) = Residual(active_dof,:);

% Concatenation of the two residues
Rf = [ R ; Ra ];

%% Forcing terms
% Should be written as if the forcing angular frequency value is 1
% i.e. the forcing period is 2*pi
Forcing = zeros(2*sys.H+1,sys.nz_tot); % DO NOT CHANGE this line.
loads = sys.parameters.loads;
dof_loc = global_to_active(loads.dof, active_dof);
Forcing(:,(dof_loc + n_var)) = loads.amplitude_cos*cos(t) + loads.amplitude_sin*sin(t);
end