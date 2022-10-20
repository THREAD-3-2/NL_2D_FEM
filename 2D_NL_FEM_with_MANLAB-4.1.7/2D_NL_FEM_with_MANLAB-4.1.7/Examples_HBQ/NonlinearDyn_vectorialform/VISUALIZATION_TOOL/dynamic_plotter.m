function dynamic_plotter(sys,U)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
set(0,'defaultTextInterpreter','latex');               %LaTeX pour le titre
set(0,'defaultAxesTickLabelInterpreter','latex');      %LaTeX pour les valeurs des axes
set(0,'defaultLegendInterpreter','latex');             %LaTeX pour les valeurs des légendes

% Defining the degrees of freedom and parameters
dof = sys.parameters.dof_info.dof;
bc_dof = sys.parameters.dof_info.bc.prescribed_dof;
active_dof = sys.parameters.dof_info.active;

number_nodes = sys.parameters.mesh.number_nodes;
number_elements = sys.parameters.mesh.number_elements;
% Le = sys.parameters.beam_params.L_beam/number_elements;

nb_period = 4;
nb_samples = 2^4;
time = linspace(0,nb_period,nb_samples);
Utime = calcperiodHBM(sys,U,[1:sys.nz/2],time);

% Initialization
q_global = zeros(3*number_nodes,1);
q_global(bc_dof) = 0; % can be set to any value
u_global = zeros(length(q_global),1);
coordinates = zeros(length(q_global),1);

for j = 1:length(time)
    q_global(active_dof) = Utime(j,:);
    for i = 1:number_elements
        nodeA = sys.parameters.mesh.connect(i,2);
        nodeB = sys.parameters.mesh.connect(i,3);

        index_global_A = 3*(nodeA-1) + [1:3];
        index_global_B = 3*(nodeB-1) + [1:3];

        index = [index_global_A index_global_B]';

        q_element = q_global(index);

        % shape functions
        N1 = 1; N2 = 0;
        N = [N1 0  0  N2 0  0;
             0  N1 0  0  N2 0;
             0  0  N1 0  0  N2;
             N2 0  0  N1 0  0;
             0  N2 0  0  N1 0;
             0  0  N2 0  0  N2];

        u_global(index) = N*q_element;
        coordinates(index_global_A') = sys.parameters.mesh.nodes(i,2:end)';
        if i == number_nodes % ring
            coordinates((3*number_nodes + [1:3])') = sys.parameters.mesh.nodes(1,2:end);
        else
            coordinates(index_global_B') = sys.parameters.mesh.nodes(i+1,2:end)';
        end
    end
%%% beam simulations
deflections = coordinates + u_global;
x_deflection = deflections(1:3:3*number_nodes);
y_deflection = deflections(2:3:3*number_nodes);

%%% ring simulations
% deflections = coordinates + [u_global; u_global(1:3)];
% x_deflection = deflections(1:3:end);
% y_deflection = deflections(2:3:end);

% Horizontal orientation
figure(12); clf; 
hold on;
plot(x_deflection,y_deflection,'b-o','LineWidth',4)
xlabel('$X$ [m]')
ylabel('$Y$[m]')
ylim([-0.5,0.5])
xlim([0,1])

% % Vertical orientation
% figure(14); clf;
% hold on;
% plot(y_deflection,x_deflection,'k-','LineWidth',4)
% xlabel('[m]')
% ylabel('[m]')
% xlim([-0.9,0.9])
% ylim([0,1])

F(j) = getframe(gcf);
end

F(j+1) = getframe(gcf);
% To create a video of figure(3):
writerObj = VideoWriter('cantilever_m3_high');
  writerObj.FrameRate = 40;
    open(writerObj);
% write the frames to the video
for i=1:length(F)
    % convert the image to a frame
    frame = F(i) ;
%     ylim([0 6])
    writeVideo(writerObj, frame);
    
end
% close the writer object
close(writerObj);
end



