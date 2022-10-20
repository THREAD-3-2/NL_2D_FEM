function [aux_expr] = auxiliary_variable_expression(obj, q_full, aux)
% expression aux auxiliary variable for the MAN equation 
% the expression must be quadratic exept for aux variables defined from a
% differential equation (see B.Cochelin MAN paper)

number_elements = obj.mesh.number_elements;
% retreive FEM variable from aux vector
up     = aux(              1:   number_elements,:);
wp     = aux(   number_elements + 1: 2*number_elements,:);
thetap = aux( 2*number_elements + 1: 3*number_elements,:);
theta  = aux( 3*number_elements + 1: 4*number_elements,:);  
c      = aux( 4*number_elements + 1: 5*number_elements,:);  
s      = aux( 5*number_elements + 1: 6*number_elements,:); 
eps    = aux( 6*number_elements + 1: 7*number_elements,:);
gam    = aux( 7*number_elements + 1: 8*number_elements,:);
Fx     = aux( 8*number_elements + 1: 9*number_elements,:); % not needed
Fy     = aux( 9*number_elements + 1:10*number_elements,:); % not needed
M      = aux(10*number_elements + 1:11*number_elements,:); % not needed
T2     = aux(11*number_elements + 1:12*number_elements,:); % not needed

for i = 1:obj.mesh.number_elements
    % Element node numbers
    nodeA  =  obj.mesh.connect(i,2);
    nodeB  =  obj.mesh.connect(i,3);
    % index of nodes dof
    index_global_A = 3*(nodeA - 1) + [1:3];
    index_global_B = 3*(nodeB - 1) + [1:3];
    index = [index_global_A index_global_B];
    % element information (length and angle w.r.t horizontal)
    [L_element, theta_element] = obj.elementOrientation(i);
    % rotation matrix to go from local frame to global frame
    rot_matrix = obj.transformationMatrix(theta_element);
    qe = rot_matrix'*q_full(index,:);    
    % gather dof from the vector qe (6x1 vector)
    u1 = qe(1,:); w1 = qe(2,:); theta1 = qe(3,:);
    u2 = qe(4,:); w2 = qe(5,:); theta2 = qe(6,:);    
    % usefull quantities (mean angle and gradient)
    aux_expr(i,:) = (u2 - u1)/L_element;       % mean axial stain
    aux_expr(i +   number_elements,:) = (w2 - w1)/L_element; % mean shear strain
    aux_expr(i +   2*number_elements,:) = (theta2 - theta1)/L_element; % mean curvature        
    aux_expr(i +   3*number_elements,:) = (theta1 + theta2)/2;   % mean theta
end

A = obj.prop.area;
I = obj.prop.inertia;
E = obj.prop.young_mod;
G = obj.prop.shear_mod;
k = obj.prop.shear_coeff_k;

N = E*A.*eps;
T = k*G*A.*gam;

% Remaining auxiliary variable expression
aux_expr( 4*number_elements + 1: 5*number_elements,:)  = cos(theta);   % c, not quadratic because differential eq             
aux_expr( 5*number_elements + 1: 6*number_elements,:)  = sin(theta);   % s, not quadratic because differential eq
aux_expr( 6*number_elements + 1: 7*number_elements,:)  = (1 + up).*c + wp.*s - 1; % eps
aux_expr( 7*number_elements + 1: 8*number_elements,:)  = wp.*c - (1 + up).*s; % gam

aux_expr( 8*number_elements + 1: 9*number_elements,:)  = N.*c - T.*s; % Fx
aux_expr( 9*number_elements + 1: 10*number_elements,:) = N.*s + T.*c; % Fy
aux_expr(10*number_elements + 1: 11*number_elements,:) = E*I*thetap;  % M
aux_expr(11*number_elements + 1: 12*number_elements,:) = N.*gam - (1 + eps).*T; % T2


end