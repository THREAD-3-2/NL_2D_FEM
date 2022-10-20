function [K_global_lin, M_global_lin, f_ext_global] = linear_matrix_assembly(parameters)

mesh = parameters.mesh;
number_elements = mesh.number_elements;
number_nodes = mesh.number_nodes;

E = parameters.beam_params.modulus;
A = parameters.beam_params.area;
G = parameters.beam_params.shear;
I = parameters.beam_params.I;
k = parameters.beam_params.k;
rho = parameters.beam_params.density;

L_element = zeros(number_elements,1);
theta_element = zeros(number_elements,1);

K_global_lin = zeros(3*number_nodes,3*number_nodes);
M_global_lin =  zeros(3*number_nodes,3*number_nodes);
f_ext_global = zeros(3*number_nodes,1);

for i = 1:number_elements
    if length(A) > 1 %%% for special case of A and I not being constant throughout the structure
        A = A(i);
        I = I(i);
    end
    nodeA  =  mesh.connect(i,2);
    nodeB  =  mesh.connect(i,3);
    
    index_global_A = 3*(nodeA - 1) + [1:3];
    index_global_B = 3*(nodeB - 1) + [1:3];

    index = [index_global_A index_global_B];
    
    [L_element(i), theta_element(i)] = elementOrientation(mesh,i);
    rot_matrix = transformationMatrix(theta_element(i));

    %%% K1 through K7 == Sénéchal's stiffness Keps + Kgam + Kkap when eps =
    %%% gam = kap = 0, c = 1, s = 0 (see Appendix of Thomas, Sénéchal, Deü 2016)
    K1 = E*A/L_element(i); 
    K3 = 0; 
    K4 = 0; 
    K2 = k*G*A/L_element(i); 
    K5 = k*G*A/2; 
    K6 = k*G*A*L_element(i)/4 + E*I/L_element(i);
    K7 = k*G*A*L_element(i)/4 - E*I/L_element(i);

    
    K_ele_lin = [ K1  K3  K4   -K1 -K3  K4;
                  K3  K2  K5   -K3 -K2  K5;
                  K4  K5  K6   -K4 -K5  K7;
                 -K1 -K3 -K4    K1  K3 -K4;
                 -K3 -K2 -K5    K3  K2 -K5;
                  K4  K5  K7   -K4 -K5  K6 ];
    
              
    M_ele_lin = (L_element(i)*rho/6)*[  2*A      0      0      A     0     0                                 
                                         0      2*A     0      0     A     0           
                                         0       0     2*I     0     0     I        
                                         A       0      0     2*A    0     0
                                         0       A      0      0    2*A    0
                                         0       0      I      0     0    2*I ];
                                     
    f_ext_elem = [0; L_element(i); 0; 0; L_element(i); 0];
    
    K_global_lin(index,index) = K_global_lin(index,index) + rot_matrix*K_ele_lin*rot_matrix';
    M_global_lin(index,index) = M_global_lin(index,index) + rot_matrix*M_ele_lin*rot_matrix';
    f_ext_global(index,:) = f_ext_global(index,:) + rot_matrix*f_ext_elem; % distributed external force
end

end