%% Plotting function

% function frfs(sys,Diagram)
% % Plot the time evolution of the mesh 
% 
set(0,'defaultTextInterpreter','latex');               %LaTeX pour le titre
set(0,'defaultAxesTickLabelInterpreter','latex');      %LaTeX pour les valeurs des axes
set(0,'defaultLegendInterpreter','latex');             %LaTeX pour les valeurs des légendes
% 
% Defining the degrees of freedom and parameters
bc_dof = sys.parameters.model.boundary.prescribed_dof;
active_dof = sys.parameters.model.boundary.active_dof;
number_nodes = sys.parameters.model.mesh.number_nodes;
% 
nb_period = 1;
nb_samples = 2^6;
time = linspace(0,nb_period,nb_samples);
% figure(30)
%%% Modify the following parameters:
mode = 1;
natural_puls = sys.parameters.natural_puls(mode);
%% Normalizing omega
for i = 1:length(Diagram)
    Diagram{1,i}.Upp(sys.neq,:) = Diagram{1,i}.Upp(sys.neq,:)/natural_puls;
end
%% FRF H of desired dof
figure(6);
a = plotbranchHBM(sys,Diagram,[sys.nz/2-1],[1],'omega'); % H1 of transverse displacement at tip
% xlabel('$\Omega/\omega_1$')
%%
omega = zeros(15*length(Diagram),1);
max_displace = zeros(15*length(Diagram),1);

for j = 1:length(Diagram)
    omega((15*(j-1) + 1):15*j) = Diagram{1,j}.Upp(sys.neq,:);
    for k = 1:15
        U = Diagram{1,j}.Upp(:,k);
        Utime = calcperiodHBM(sys,U,[1:sys.nz/2],time);
        maxTransverse = max(abs(Utime(:,2:3:end)));
        [rowOfMax, colOfMax] = find(abs(Utime) == max(maxTransverse));
        
        max_displace(15*(j-1) + k) = (max(Utime(:,colOfMax(1))) - min(Utime(:,colOfMax(1))))/2;
%         max_displace(15*(j-1) + k) = max(U(1:(sys.nz/2)*(2*sys.H + 1)));
%         max_displace(15*(i-1) + j) = max(U(2:3:(sys.nz/2)*(2*sys.H + 1)));
% need to reconstruct
    end
end
figure(7);
plot(omega,max_displace)
xlim([0.85,1.15])
% end



