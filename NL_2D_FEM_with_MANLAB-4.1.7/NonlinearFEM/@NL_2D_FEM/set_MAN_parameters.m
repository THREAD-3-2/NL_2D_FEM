function [nz, nz_aux, parameters] = set_MAN_parameters(obj, H, type, model, angfreq)
% creates the 'parameters' structure for the MAN computation
% general parameters
parameters.type = type;
parameters.model = model;
parameters.H = H;
parameters.angfreq = angfreq; % 'omega' if the frequency is a parameter, 
                              %  value  if the system is forced at a fixed angular frequency
% Parameters relative to the number of equation
neq = length(obj.boundary.active_dof); %% number of variables: [u,w,theta] per node
neq_aux = 12*obj.mesh.number_elements; %% number of auxiliary variables: in this case, 12 per element
nz = 2*neq;         % number of main equations of the system of the differential-algebraic system (DAE)
nz_aux = neq_aux;   % number of auxiliary equations of the system of the DAE

end
