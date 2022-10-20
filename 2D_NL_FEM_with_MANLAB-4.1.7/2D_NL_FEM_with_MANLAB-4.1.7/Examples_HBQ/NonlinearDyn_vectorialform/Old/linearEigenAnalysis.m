function [parameters, K_global_lin, M_global_lin, K, M, C, omega, V, omega_free, Phi] = linearEigenAnalysis(parameters)
% K_glob, M_glob full matrices
% M, K, C matrices taking into acount BC (reduced size)
% Phi, omega_free : free modes
% V, omega : modes taking into account BC

%%%
[K_global_lin, M_global_lin, f_ext_global] = linear_matrix_assembly(parameters);

% free linear modes
[Phi,D] = eig(K_global_lin, M_global_lin);
[omega_free, idx] = sort(sqrt(diag(D)), 'ascend');
Phi = Phi(:,idx);

% modal analysis with boundary conditions
prescribed_dof = parameters.dof_info.bc.prescribed_dof;
M = M_global_lin;
K = K_global_lin;
M(prescribed_dof,:) = [ ];  % remove rows with fixed dof
M(:,prescribed_dof) = [ ];  % remove columns with fixed dof
K(prescribed_dof,:) = [ ];
K(:,prescribed_dof) = [ ];
f_ext_global(prescribed_dof) = [ ];


[V, D] = eig(K, M); 
[omega, idx] = sort(sqrt(diag(D)), 'ascend');
V = V(:,idx);

% linear damping matrix
if strcmp(parameters.beam_params.damping.type,'Rayleigh')
%     C = parameters.beam_params.damping.alpha * M + parameters.beam_params.damping.beta * K;
   alpha = 2 * parameters.beam_params.damping.xi0 * omega(parameters.mode);
   beta = 0;
   C = alpha * M + beta * K;
elseif strcmp(parameters.beam_params.damping.type,'Modal')
%     C = M * V * parameters.beam_params.damping.Xi * transpose(V) * M;
    parameters.beam_params.damping.Xi(parameters.mode,parameters.mode) = parameters.beam_params.damping.Xi_mode; % reduce damping at selected mode
%     parameters.beam_params.damping.Xi(2:11,2:11) = parameters.beam_params.damping.Xi_increased*eye(10); % increase damping on first 10 transverse modes only
%     parameters.beam_params.damping.Xi = parameters.beam_params.damping.Xi./(2*omega); % checking that Modal = Raleigh
    C = M * V * 2 * parameters.beam_params.damping.Xi * diag(omega) * transpose(V) * M;
else
    disp('error: undefined damping type')
    return
end

% recast to full size, with zeros at BC dof position
C_global_lin = zeros(size(M_global_lin));
C_global_lin(parameters.dof_info.active, parameters.dof_info.active) = C;

% put the matrices inside the structure
parameters.M_glob = M_global_lin;
parameters.K_glob = K_global_lin;
parameters.C_glob = C_global_lin;
parameters.f_ext_glob = f_ext_global;
end

