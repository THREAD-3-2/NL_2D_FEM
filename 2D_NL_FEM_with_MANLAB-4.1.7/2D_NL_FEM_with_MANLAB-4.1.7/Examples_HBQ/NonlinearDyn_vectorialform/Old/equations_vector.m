function [Rf,dRa,Forcing] = equations_vector(sys, t, zf, dzf, d2zf)
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
rho = sys.parameters.beam_params.density;
%% Variables : 
u      = zf(1:n_var,:);           du     = dzf(1:n_var,:);             % Main variable: u = position
v      = zf(n_var + 1:2*n_var,:); dv     = dzf(n_var + 1:2*n_var,:);   % Main variable: v = velocity
aux    = zf(2*n_var + 1:end-1,:); daux   = dzf(2*n_var + 1:end,:);     % Auxiliary variables
lambda = zf(end,:);                                                  % Continuation parameter

% decomposition of aux:
up     = aux(              1:   number_elements,:);
wp     = aux(   number_elements + 1: 2*number_elements,:);
thetap = aux( 2*number_elements + 1: 3*number_elements,:);
theta  = aux( 3*number_elements + 1: 4*number_elements,:);  dtheta = daux(3*number_elements + 1:4*number_elements,:);
c      = aux( 4*number_elements + 1: 5*number_elements,:);      dc = daux(4*number_elements + 1:5*number_elements,:);
s      = aux( 5*number_elements + 1: 6*number_elements,:);      ds = daux(5*number_elements + 1:6*number_elements,:);
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
f_ext_global = sys.parameters.f_ext_glob;

for i = 1:number_elements
    nodeA = sys.parameters.mesh.connect(i,2);
    nodeB = sys.parameters.mesh.connect(i,3);
    
    index_global_A = 3*(nodeA - 1) + (1:3);
    index_global_B = 3*(nodeB - 1) + (1:3);
    
    index = [index_global_A index_global_B]';
    
    [L_element(i), theta_element(i)] = elementOrientation(sys.parameters.mesh,i);
    rot_matrix = transformationMatrix(theta_element(i));
    q_element = rot_matrix'*q_global(index,:); % rotating q_element into local frame
    
    f_int_elem = [-Fx(i,:);
        -Fy(i,:);
        -M(i,:) + T2(i,:)*L_element(i)/2; 
        Fx(i,:);
        Fy(i,:);
        M(i,:) + T2(i,:)*L_element(i)/2];
    f_int_global(index,:) = f_int_global(index,:) + rot_matrix*f_int_elem;
    
%     if strcmp(parameters.loads.type, 'Distributed force')
%         f_ext_elem = [0; L_element(i)/2; 0; 0; L_element(i); 0];
%         f_ext_global(index,:) = f_ext_global(index,:) + rot_matrix*f_ext_elem; % distributed external force
%     end

    Ra(i,:)                 = up(i,:)     - ((q_element(4,:) - q_element(1,:))/L_element(i));   
    Ra(i +   number_elements,:) = wp(i,:)     - ((q_element(5,:) - q_element(2,:))/L_element(i));
    Ra(i + 2*number_elements,:) = thetap(i,:) - ((q_element(6,:) - q_element(3,:))/L_element(i));
    Ra(i + 3*number_elements,:) = theta(i,:)  - ((q_element(6,:) + q_element(3,:))/2);
    
end
N = E*A.*eps;
T = k*G*A.*gam;

% Remaining auxiliary variable residuals
Ra( 4*number_elements + 1: 5*number_elements,:)  = c - cos(theta);                 dRa(4*number_elements + 1: 5*number_elements,:) = dc + s.*dtheta;
Ra( 5*number_elements + 1: 6*number_elements,:)  = s - sin(theta);                 dRa(5*number_elements + 1: 6*number_elements,:) = ds - c.*dtheta;
Ra( 6*number_elements + 1: 7*number_elements,:)  = eps - (1 + up).*c - wp.*s + 1;
Ra( 7*number_elements + 1: 8*number_elements,:)  = gam - wp.*c + (1 + up).*s;

Ra( 8*number_elements + 1: 9*number_elements,:)  = Fx - N.*c + T.*s;
Ra( 9*number_elements + 1: 10*number_elements,:) = Fy - N.*s - T.*c;
Ra(10*number_elements + 1: 11*number_elements,:) = M - E*I*thetap;
Ra(11*number_elements + 1: 12*number_elements,:) = T2 - N.*gam + (1 + eps).*T;

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
if strcmp(loads.type, 'Point force') % point force at specific nodes
    Forcing(:,(dof_loc + n_var)) = loads.amplitude_cos*cos(t) + loads.amplitude_sin*sin(t);
elseif strcmp(loads.type, 'Distributed force') % distributed force due to imposed acceleration
%     Forcing(:,(dof_loc + n_var)) = rho * A * repmat(f_ext_global',size(Forcing(:,(dof_loc + n_var)),1),1) * (loads.amplitude_cos*cos(t) + loads.amplitude_sin*sin(t));
    for i = 1:length(dof_loc)
        Forcing(:,(dof_loc(i) + n_var)) = rho * A * f_ext_global(i) * (loads.amplitude_cos*cos(t) + loads.amplitude_sin*sin(t));
    end
end
end