function visualize_point(sys,fig,Defs)
%
%   
set(0,'defaultTextInterpreter','latex');               %LaTeX pour le titre
set(0,'defaultAxesTickLabelInterpreter','latex');      %LaTeX pour les valeurs des axes
set(0,'defaultLegendInterpreter','latex');             %LaTeX pour les valeurs des légendes

k = [2 4]; % identifying location of subplots in visualization plot
Names = fieldnames(Defs);

% Defining the degrees of freedom and parameters
bc_dof = sys.parameters.model.boundary.prescribed_dof;
active_dof = sys.parameters.model.boundary.active_dof;
number_nodes = sys.parameters.model.mesh.number_nodes;

nb_period = 1;
nb_samples = 6;
time = linspace(0,nb_period,nb_samples);
for jj = 1:numel(Names)
    U = Defs.(Names{jj,:});
    Utime = calcperiodHBM(sys,U,[1:sys.nz/2],time);

    % Initialization
    q_global = zeros(3*number_nodes,1);
    q_global(bc_dof) = 0; % can be set to any value

    figure(fig)
    subplot(2,2,k(jj)) % k indicates which subplot to plot in

    if isfield(sys.parameters.model.visu.axis_info,'axis')
        axis(sys.parameters.dof_info.axis)
        axis square
    else
        axis([-1.5, 1.5, -1.5, 1.5])
        axis square
    end
    for j = 1:length(time)
        q_global(active_dof) = Utime(j,:);
    %cla
        plot_deformed_mesh(sys.parameters.mesh,q_global,fig,'-b') % plot deformed mesh at time tj
    end
    plot_deformed_mesh(sys.parameters.mesh,0*q_global,fig,'--k') % plot deformed mesh at time tj
    xlabel('$X$ [m]')
    ylabel('$Y$ [m]')

    % Plot point on FRF
    mode = sys.parameters.mode; % mode being visualized
    natural_puls = sys.parameters.natural_puls; % natural frequencies
    harm = 1; % amplitude of this harmonic will be plotted
    H = sys.parameters.H; % number of harmonics used to compute the solution
    obs_dof = sys.parameters.dof_info.obs_dof; % degree of freedom to be observed (discarding BC dof)

    % Finding locations in high amplitude and low amplitude 'Point's
    I0 = (obs_dof - 1)*(2*H + 1) + 1; % constant term in harmonic expansion
    Icos = (obs_dof - 1)*(2*H + 1) + 1 + harm; % coefficient in front cosine of harmonic to be plotted
    Isin = (obs_dof - 1)*(2*H + 1) + 1 + H + harm; % coefficient in front sine of harmonic to be plotted

    harm_amp_high = sqrt(U(Icos)^2 + U(Isin)^2);
    omega_high = U(sys.neq)/natural_puls(mode);

    figure(fig)
    subplot(2,2,[1 3])
    % plot(omega_high,harm_amp_high,'ko','LineWidth',2,'DisplayName','High Amp. Point')
    plot(omega_high,harm_amp_high,'o','LineWidth',2)
end
end

