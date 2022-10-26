function [strain, stress] = strains_and_stress_at_gauss_point(obj, q_full)
% computes stress and strains

% element material properties
A = obj.prop.area;
I = obj.prop.inertia;
E = obj.prop.young_mod;
G = obj.prop.shear_mod;
k = obj.prop.shear_coeff_k;
% initialise lists
Ncol = size(q_full,2);
theta = zeros(obj.mesh.number_elements,Ncol);
thetap = zeros(obj.mesh.number_elements,Ncol);
up = zeros(obj.mesh.number_elements,Ncol);
wp = zeros(obj.mesh.number_elements,Ncol);
% loop over elements
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
    theta(i,:) = (theta1 + theta2)/2;   % mean theta
    up(i,:) = (u2 - u1)/L_element;       % mean axial stain
    wp(i,:) = (w2 - w1)/L_element; % mean shear strain
    thetap(i,:) = (theta2 - theta1)/L_element; % mean curvature        
end
% fill outputs
strain.meanrot = theta;
strain.grad.up = up;
strain.grad.wp = wp;
strain.grad.thetap = thetap;
% compute remaining auxiliary variable 
strain.c = cos(theta);                 
strain.s = sin(theta);
% strain
strain.eps = (1 + up).*strain.c + wp.*strain.s - 1;   
strain.gam = wp.*strain.c - (1 + up).*strain.s;
strain.chi = thetap;
% stress
stress.sigma = E*strain.eps ;
stress.tau =   G*strain.gam;
% internal forces (local frame)
stress.N = E*A.*strain.eps;
stress.T = k*G*A.*strain.gam;
stress.M = E*I*thetap;                 
% internal forces (global frame)
stress.Fx = stress.N.*strain.c - stress.T.*strain.s;               
stress.Fy = stress.N.*strain.s + stress.T.*strain.c;
stress.T2 = stress.N.*strain.gam - (1 + strain.eps).*stress.T;
end