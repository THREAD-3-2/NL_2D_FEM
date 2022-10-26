function plot_deformed_mesh_with_stress(obj, q_global, stress, fig, marker)
% plot the deformed configuration contained in q_glob
% q_glob is of full size and contain 0's at the position of prescribed dof

figure(fig);
hold on;

N  = stress.N;
T  = stress.T;
M  = stress.M;

% Sress boundary
Nmax = max(max(N),1); Nmin = min(0,min(N));
Tmax = max(max(T),1); Tmin = min(0,min(T));
Mmax= max(max(M),1);  Mmin= min(0,min(M));

fprintf(1,'\nNmax : %2.2e , Nmin : %2.2e\n', Nmax, Nmin)
fprintf(1,'Tmax : %2.2e , Tmin : %2.2e\n', Tmax, Tmin)
fprintf(1,'Mmax : %2.2e , Mmin : %2.2e\n \n', Mmax, Mmin)

% mesh boundary
maxX = max(obj.mesh.nodes(:,2));
maxY = max(obj.mesh.nodes(:,3));
minX = min(obj.mesh.nodes(:,2));
minY = min(obj.mesh.nodes(:,3));

% color map
 cmap = colormap('parula'); % take your pick (doc colormap)
% cmap = colormap('jet'); % take your pick (doc colormap)
%
figN = figure(666); figN.Name='Normal stress'; 
figT = figure(667); figT.Name='Transverse stress'; 
figM = figure(668); figM.Name='Bending torque'; 
fig = {figN, figT, figM};

for i = 1:obj.mesh.number_elements
    nodeA = obj.mesh.connect(i,2);  % first node of element
    nodeB = obj.mesh.connect(i,3);  % second node of element
    
    index_global_A = 3*(nodeA-1) + [1:3]; % index of dof of node A
    index_global_B = 3*(nodeB-1) + [1:3]; % index of dof of node B
    
    index = [index_global_A index_global_B]';
    
    q_e = q_global(index); % element displacement vector
    
    P1 = obj.mesh.nodes(nodeA, 2:3) + q_e(1:2)'; % position of node A after deformation
    P2 = obj.mesh.nodes(nodeB, 2:3) + q_e(4:5)'; % position of node B after deformation
    
    for pp=1:3
        figure(fig{pp}); hold on
        h = plot([P1(1), P2(1)],[P1(2) P2(2)], marker, 'linewidth', 1); % plot the current element (deformed)
        plot(P1(1), P1(2),'.k', 'markersize', 8) % plot the nodes (in deformed position)
        plot(P2(1), P2(2),'.k', 'markersize', 8)
        
        if pp==1
            N(i);
            cint = interp1(linspace(Nmin,Nmax,length(cmap)),cmap,N(i)); % map color to y values
        elseif pp==2
            T(i);
            cint = interp1(linspace(Tmin,Tmax,length(cmap)),cmap,T(i)); % map color to y values
        elseif pp==3
            M(i);
            cint = interp1(linspace(Mmin,Mmax,length(cmap)),cmap,M(i)); % map color to y values
        end
        cint = uint8(cint'*255); % need a 4xN uint8 array
        cint(4,:) = 255; % last column is transparency
        set(h,'Color',cint)
    end
end
    drawnow
    
    for pp=1:3
        figure(fig{pp}); c = colorbar;
        c.TickLabelsMode = 'manual';
        c.TicksMode = 'manual';

        if pp==1
        c.Label.String = 'Normal stress';
        c.Ticks = linspace(Nmin, Nmax, 5);
        c.TickLabels = linspace(Nmin, Nmax, 5);
        c.Limits = [Nmin, Nmax];

        elseif pp==2
        c.Label.String = 'Transverse stress';
        c.Ticks = linspace(Tmin, Tmax, 5);
        c.TickLabels = linspace(Tmin, Tmax, 5);
        c.Limits = [Tmin, Tmax];
        elseif pp==3
        c.Label.String = 'Bending moment';
        c.Ticks = linspace(Mmin, Mmax, 5);
        c.TickLabels = linspace(Mmin, Mmax, 5);
        c.Limits = [Mmin, Mmax];
        end

%        keyboard
    end
end