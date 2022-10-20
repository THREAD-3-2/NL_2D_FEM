function [mesh, mother_geometric_line] = mesh_from_truss(obj, geometricPoint, geometricElement, geometricDiscretisation)

coord = geometricPoint; % add geometric point to the FEM nodes
Np0 = size(coord,1);

Ne0 = 0;      % initialiase number of element to zero
element = []; % table of connectivity to be constructed
mother_geometric_line = []; % list of line (geometric element) number to witch the (FEM) elements belongs
tol = 1e-4;   % tolerence for point detection

for i=1:size(geometricElement,1)
    A = geometricElement(i,2); % index of first point
    B = geometricElement(i,3); % index of second point
    P1 = geometricPoint(A,2:3); % coord of first point
    P2 = geometricPoint(B,2:3); % coord of second point
    
    % search if the geometric points are already in the list of nodes
    numP1 = find(sqrt((P1(1)-coord(:,2)).^2+(P1(2)-coord(:,3)).^2)<tol, 1, 'first');
    numP2 = find(sqrt((P2(1)-coord(:,2)).^2+(P2(2)-coord(:,3)).^2)<tol, 1, 'first');
    
    % nuber of (discretized) element for the curent geometric element
    Ne = geometricDiscretisation(i);
    % update element list with the element of the current line
    [coord, element, Np0, Ne0] = mesh_line(P1, P2, Ne, Ne0, Np0, coord, element, numP1, numP2);
    
    mother_geometric_line = [mother_geometric_line, repmat(i,1,Ne)];
end
mesh.nodes = coord;
mesh.connect = element;
number_nodes = size(coord,1);
mesh.number_nodes = number_nodes;
mesh.number_elements = size(element,1); 
% mesh.mother_geometric_line = mother_geometric_line;
end

function [coord, element, Np0, Ne0] = mesh_line(P1, P2, Ne, Ne0, Np0, coord, element, numP1, numP2)
Np00 = Np0;
% parametrisation du segment
% P = (1-lambda) P1 + lambda P2  avec lambda = l/L in [0, 1]
lambda = linspace(0,1,Ne+1);

if Ne>1
% Deal with coord list
for i = 2:length(lambda)-1 % discard boundary point (already in the coord list from initialisation)
    Np0 = Np0+1;           % add one point
    P = (1-lambda(i))*P1 + lambda(i)*P2;
    coord_p = [ Np0 P(1) P(2) 0];
    coord = [coord; coord_p]; % update the coordinates list
end
end

% Deal with element (i.e. connectivity table) list
if Ne==1
Ne0 = Ne0+1; % add first element of the line
element_c = [Ne0 numP1 numP2];
element = [element; element_c];
else
Ne0 = Ne0+1; % add first element of the line
element_c = [Ne0 numP1 Np00+1];
element = [element; element_c];
end
if Ne>=3
for i = 1:Ne-2
    Ne0 = Ne0+1; % add one element
    element_c = [Ne0 Np00+1 Np00+2];
    Np00 = Np00+1;
    element = [element; element_c];
end
end
if Ne>=2
Ne0 = Ne0+1; % add last element of the line
element_c = [Ne0 Np00+1 numP2];
element = [element; element_c];
end

end
