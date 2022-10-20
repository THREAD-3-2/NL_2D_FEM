function [Xhat, AmpH1, PhaseH1] = linear_solution(M, C, K, Omega, parameters, obs_dof)
% Real(Xhat) = XcosH1
% Imag(Xhat) = XsinH1
% obs_dof = observed dof (discarding BC dof)
ndof = size(M,1);
Nw = length(Omega);

Xhat = zeros(ndof,Nw);
AmpH1 = zeros(ndof,Nw);
PhaseH1 = zeros(ndof,Nw);

% create force vector 
loads = parameters.loads;
dof_loc = global_to_active(loads.dof, parameters.dof_info.active);
Fhat = zeros(ndof,1);
if strcmp(loads.type, 'Point force') % point force at specific nodes
    Fhat(dof_loc) = loads.amplitude_cos + 1i*loads.amplitude_sin;
elseif strcmp(loads.type, 'Distributed force') % distributed force due to imposed acceleration
    Fhat(dof_loc) = parameters.beam_params.density * parameters.beam_params.area * parameters.f_ext_glob *(loads.amplitude_cos + 1i*loads.amplitude_sin);
end

for i=1:Nw
   Xhat(:,i) = (K-Omega(i)^2*M+1i*C*Omega(i))\Fhat;
   AmpH1(:,i) = abs(Xhat(:,i));
   PhaseH1(:,i) = angle(Xhat(:,i));
end

figure(666)
subplot(2,1,1)
plot(Omega, AmpH1(obs_dof,:))
grid on
xlabel('omega')
ylabel('Amplitude')
subplot(2,1,2)
plot(Omega, unwrap(PhaseH1(obs_dof,:)))
grid on
xlabel('omega')
ylabel('Phase')

figure(5)
plot(Omega, AmpH1(obs_dof,:))
grid on
xlabel('omega')
ylabel('Amplitude')
end