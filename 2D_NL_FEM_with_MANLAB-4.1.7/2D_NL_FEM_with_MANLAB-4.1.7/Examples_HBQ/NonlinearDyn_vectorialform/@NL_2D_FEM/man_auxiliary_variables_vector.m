function [aux] = man_auxiliary_variables_vector(obj, q_full)
[strain, stress] = strains_and_stress_at_gauss_point(obj, q_full);
number_elements = obj.mesh.number_elements;
% decomposition of aux:
aux(              1:   number_elements,:) = strain.grad.up;
aux(   number_elements + 1: 2*number_elements,:) = strain.grad.wp;
aux( 2*number_elements + 1: 3*number_elements,:) = strain.grad.thetap;
aux( 3*number_elements + 1: 4*number_elements,:) = strain.meanrot;  
aux( 4*number_elements + 1: 5*number_elements,:) = strain.c;      
aux( 5*number_elements + 1: 6*number_elements,:) = strain.s;      
aux( 6*number_elements + 1: 7*number_elements,:) = strain.eps;
aux( 7*number_elements + 1: 8*number_elements,:) = strain.gam;
aux( 8*number_elements + 1: 9*number_elements,:) = stress.Fx;
aux( 9*number_elements + 1:10*number_elements,:) = stress.Fy;
aux(10*number_elements + 1:11*number_elements,:) = stress.M;
aux(11*number_elements + 1:12*number_elements,:) = stress.T2;
end