function [Z0,ModalDisp] = nnmInitialization(sys,nz,nz_aux,H,omega0,modeShape,epsilon)
% nnmInitialization.m outputs an initial Z0 matrix of variables not far
% from the solution, making it more likely that MANLAB will find and trace
% the nonlinear normal modes of the system

ModalDisp = zeros(length(sys.parameters.dof_info.dof),1);
ModalDisp(sys.parameters.dof_info.active) = modeShape*epsilon;

number_elements = sys.parameters.mesh.number_elements;

E = sys.parameters.beam_params.modulus;
A = sys.parameters.beam_params.area;
I = sys.parameters.beam_params.I;
G = sys.parameters.beam_params.shear;

% compute initial gradient (in local coordinates, uses the rotation matrix)
for i = 1:number_elements
    nodeA = sys.parameters.mesh.connect(i,2);
    nodeB = sys.parameters.mesh.connect(i,3);    
    
    index_global_A = 3*(nodeA - 1) + (1:3);
    index_global_B = 3*(nodeB - 1) + (1:3);
    index = [index_global_A index_global_B]';
    
    [L_element(i), theta_element(i)] = elementOrientation(sys.parameters.mesh,i);
    rot_matrix = transformationMatrix(theta_element(i));
    q_element = rot_matrix'*ModalDisp(index); % rotating q_element into local frame
    
    up(i)    = ((q_element(4,:) - q_element(1,:))/L_element(i));   
    wp(i)    = ((q_element(5,:) - q_element(2,:))/L_element(i));
    thetap(i)= ((q_element(6,:) - q_element(3,:))/L_element(i));
    theta(i) = ((q_element(6,:) + q_element(3,:))/2);    
end

% initialization of the auxiliary variables
C   =  cos(theta);
S   =  sin(theta);
eps =  (up);
gam =  (wp - theta);
N   =  E*A*eps;
T   =  G*A*gam;
Fx  =  (N);
Fy  =  (T);
M   =  E*I*thetap;
T2  =  -T;

Z0 = 1e-7*rand(2*H+1, sys.nz_tot);
% The matrix Z0 contains in column the Fourier development of all the
% variables of your system :
% Z0 = [ u , v ]
% where u = [ u_0 ; u_c1 ; u_c2 ; ... ; u_cH ; u_s1 ; u_s2 ; ... ; u_sH ];
% u_0 is the constant Fourier coefficient, u_c1 the first cosine Fourier
% coefficient, u_s1 the first sine Fourier coefficient, etc...
% Z0(1,1)     = u0;
Z0(2,1:nz/2) = ModalDisp(sys.parameters.dof_info.active)'; % initial position
Z0(2+H,nz/2+1:nz) = -omega0*ModalDisp(sys.parameters.dof_info.active)'; % initial position

Z0(2,nz + (1:number_elements)) = up; % equations_vector.m
Z0(2,nz + (number_elements + 1:2*number_elements))     = wp;     % equations_vector.m
Z0(2,nz + (2*number_elements + 1:3*number_elements))   = thetap; % equations_vector.m
Z0(2,nz + (3*number_elements + 1:4*number_elements))   = theta;  % equations_vector.m

Z0(1,nz + (4*number_elements + 1:5*number_elements))   = C;      % theta at gauss point
Z0(2,nz + (4*number_elements + 1:5*number_elements))   = 0;      % theta at gauss point

Z0(2,nz + (5*number_elements + 1:6*number_elements))   = S;      % theta at gauss point
Z0(2,nz + (6*number_elements + 1:7*number_elements))   = eps;    % theta at gauss point
Z0(2,nz + (7*number_elements + 1:8*number_elements))   = gam;    % theta at gauss point
Z0(2,nz + (8*number_elements + 1:9*number_elements))   = Fx;     % theta at gauss point
Z0(2,nz + (9*number_elements + 1:10*number_elements))  = Fy;     % theta at gauss point
Z0(2,nz + (10*number_elements + 1:11*number_elements)) = M;      % theta at gauss point
Z0(2,nz + (11*number_elements + 1:12*number_elements)) = T2;     % theta at gauss point

end

