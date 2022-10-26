function Kte = elementary_tangent_stiffness_matrix(obj, L_element, qe)
% elementary tangent matrix of the geometrically exact element (see
% A.Senechal paper). qe is the element dof vector of size 6x1

% element material properties
A = obj.prop.area;
I = obj.prop.inertia;
E = obj.prop.young_mod;
G = obj.prop.shear_mod;
k = obj.prop.shear_coeff_k;

% gather dof from the vector qe (6x1 vector)
u1 = qe(1); w1 = qe(2); theta1 = qe(3);
u2 = qe(4); w2 = qe(5); theta2 = qe(6);

% usefullquantities
theta_bar = (theta1 + theta2)/2;   % mean theta
e_bar = (1+(u2 - u1)/L_element)*cos(theta_bar) + ((w2-w1)/L_element)*sin(theta_bar) -1;       % mean axial stain
gam_bar = (w2 - w1)/L_element*cos(theta_bar) - (1+(u2 - u1)/L_element)*sin(theta_bar); % mean shear strain
kap_bar = (theta2 - theta1)/L_element; % mean curvature

c_bar = cos(theta_bar);
s_bar = sin(theta_bar);

K1 = c_bar^2;
K2 = s_bar^2;
K3 = c_bar*s_bar;
K4 = L_element/2*( e_bar*s_bar - gam_bar*c_bar);
K5 = -L_element/2*(e_bar*c_bar + gam_bar*s_bar);
K6 = L_element^2/4*(gam_bar^2 - e_bar*(e_bar+1));
K7 = L_element/2*(gam_bar*c_bar - (e_bar+1)*s_bar);
K8 = L_element/2*(gam_bar*s_bar + (e_bar + 1)*c_bar);
K9 = L_element^2/4*((e_bar + 1)^2-gam_bar^2);

Ke_element = E*A/L_element*[ K1  K3  K4 -K1 -K3  K4;
    K3  K2  K5 -K3 -K2  K5;
    K4  K5  K6 -K4 -K5  K6;
    -K1 -K3 -K4  K1  K3 -K4;
    -K3 -K2 -K5  K3  K2 -K5;
    K4  K5  K6 -K4 -K5  K6 ]; % axial part

Kgam_element = k*G*A/L_element*[ K2 -K3  K7 -K2  K3  K7;
    -K3  K1  K8  K3 -K1  K8;
    K7  K8  K9 -K7 -K8  K9;
    -K2  K3 -K7  K2 -K3 -K7;
    K3 -K1 -K8 -K3  K1 -K8;
    K7  K8  K9 -K7 -K8  K9 ]; % transverse part

Kkap_element = E*I/L_element*[ 0 0  0 0 0  0;
    0 0  0 0 0  0;
    0 0  1 0 0 -1;
    0 0  0 0 0  0;
    0 0  0 0 0  0;
    0 0 -1 0 0  1 ]; % rotation part

Kte = Ke_element + Kgam_element + Kkap_element;


%  K1 = E*A/L_element; 
%     K3 = 0; 
%     K4 = 0; 
%     K2 = k*G*A/L_element; 
%     K5 = k*G*A/2; 
%     K6 = k*G*A*L_element/4 + E*I/L_element;
%     K7 = k*G*A*L_element/4 - E*I/L_element;
% 
%     
%     Kte = [ K1  K3  K4   -K1 -K3  K4;
%                   K3  K2  K5   -K3 -K2  K5;
%                   K4  K5  K6   -K4 -K5  K7;
%                  -K1 -K3 -K4    K1  K3 -K4;
%                  -K3 -K2 -K5    K3  K2 -K5;
%                   K4  K5  K7   -K4 -K5  K6 ];

end
