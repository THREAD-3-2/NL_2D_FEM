function Fe = elementary_internal_force_vector(obj, qe, L_element)
% vector of internal forces for the geometrically exact 2D FE 
A = obj.prop.area;
I = obj.prop.inertia;
E = obj.prop.young_mod;
G = obj.prop.shear_mod;
k = obj.prop.shear_coeff_k;

% gather dof from the vector qe (6x1 vector)
u1 = qe(1,:); w1 = qe(2,:); theta1 = qe(3,:);
u2 = qe(4,:); w2 = qe(5,:); theta2 = qe(6,:);

% usefull quantities
theta = (theta1 + theta2)/2;   % mean theta
up = (u2 - u1)/L_element;       % mean axial stain
wp = (w2 - w1)/L_element; % mean shear strain
thetap = (theta2 - theta1)/L_element; % mean curvature

% Remaining auxiliary variable residuals
c = cos(theta);                 
s = sin(theta);
eps = (1 + up).*c + wp.*s - 1;   
gam = wp.*c - (1 + up).*s;

N = E*A.*eps;
T = k*G*A.*gam;

Fx = N.*c - T.*s;               
Fy = N.*s + T.*c;
M = E*I*thetap;                 
T2 = N.*gam - (1 + eps).*T;

Fe = [-Fx;
      -Fy;
        -M + T2*L_element/2;
       Fx;
       Fy;
        M + T2*L_element/2];

end
