function harmonicDecomposition(sys,U)
% Fourier coefficients of selected variables at point U
% Currently plots harmonic decomposition of auxiliary variables c, s, eps
% and gam at selected node (List of aux variables: [up wp thetap theta c s
% eps gam Fx Fy M T2]
number_elements = sys.parameters.mesh.number_elements;
obs_dof = sys.parameters.dof_info.obs_dof; % observed dof in FRF/NNM calculation
obs_node = ceil(obs_dof/3); % node at which harmonic decomposition is plotted
H = sys.parameters.H; % number of harmonics kept in Fourier decomposition
vars = sys.nz + (4:7)*number_elements + obs_node; % c, s, eps, gam aux variables at obs_node

amps = zeros(4*(H + 1)); % initializing vector of coefficient amplitudes
harmonics = 0:H;
for jj = 1:length(vars)
    I0 = (vars(jj) - 1)*(2*H + 1) + 1 + 2;
    amps((H + 1)*(jj - 1) + 1) = U(I0);
    
    Icos = (vars(jj) - 1)*(2*H + 1) + 1 + (1:H); 
    Isin = (vars(jj) - 1)*(2*H + 1) + 1 + H + (1:H);
    amps(((H + 1)*(jj - 1) + 2):((H + 1)*jj)) = sqrt(U(Icos).^2 + U(Isin).^2);
end
c_amp = amps(1:(H+1));
s_amp = amps((H+1)+1:2*(H+1));
eps_amp = amps(2*(H+1)+1:3*(H+1));
gam_amp = amps(3*(H+1)+1:4*(H+1));

% plot_amp = [c_amp; s_amp; eps_amp; gam_amp]';

figure
bar(harmonics,[c_amp; s_amp]')
xlabel('Harmonic')
ylabel('Amplitude')
legend('$c = \cos(\theta)$','$s = \sin(\theta)$')

figure
bar(harmonics,[eps_amp; gam_amp]')
xlabel('Harmonic')
ylabel('Amplitude')
legend('$e$','$\gamma$')

figure
bar(harmonics,[eps_amp.*gam_amp]')
xlabel('Harmonic')
ylabel('Amplitude')
legend('$e*\gamma$')
end

