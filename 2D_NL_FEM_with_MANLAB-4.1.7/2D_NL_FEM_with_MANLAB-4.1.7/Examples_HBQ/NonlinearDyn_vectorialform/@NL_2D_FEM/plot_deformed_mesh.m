function plot_deformed_mesh(obj, q_global, fig, marker)
% plot the deformed configuration contained in q_glob
% q_glob is of full size and contain 0's at the position of prescribed dof

figure(fig);
hold on;

maxX = max(obj.mesh.nodes(:,2));
maxY = max(obj.mesh.nodes(:,3));
minX = min(obj.mesh.nodes(:,2));
minY = min(obj.mesh.nodes(:,3));

for i = 1:obj.mesh.number_elements
    nodeA = obj.mesh.connect(i,2);  % first node of element
    nodeB = obj.mesh.connect(i,3);  % second node of element
    
    index_global_A = 3*(nodeA-1) + [1:3]; % index of dof of node A
    index_global_B = 3*(nodeB-1) + [1:3]; % index of dof of node B
    
    index = [index_global_A index_global_B]';
    
    q_e = q_global(index); % element displacement vector
    
    P1 = obj.mesh.nodes(nodeA, 2:3) + q_e(1:2)'; % position of node A after deformation
    P2 = obj.mesh.nodes(nodeB, 2:3) + q_e(4:5)'; % position of node B after deformation
    
    plot([P1(1), P2(1)],[P1(2) P2(2)], marker, 'linewidth', 1) % plot the current element (deformed)
    plot(P1(1), P1(2),'.k', 'markersize', 8) % plot the nodes (in deformed position)
    plot(P2(1), P2(2),'.k', 'markersize', 8)
%    keyboard
end


