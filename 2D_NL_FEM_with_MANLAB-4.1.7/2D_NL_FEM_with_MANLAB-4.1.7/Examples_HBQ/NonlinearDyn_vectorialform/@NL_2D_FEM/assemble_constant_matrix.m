function [A_full] = assemble_constant_matrix(obj, type, varargin)
% assemble a matrix depending on the provided type
% type: 'mass', 'stiffness_at_origin', 'stiffness_at_qs'
% varargin depending on type
% mass or stiffness_at_origin : varargin empty
% stiffness_at_qs: varargin{1} = qs (given configuration vector)
%

% deal with variable arguments
if strcmp(type,'stiffness_at_qs')
    qs = varargin{1} ;
end
% initialize matrix
A_full = zeros(3*obj.mesh.number_nodes, 3*obj.mesh.number_nodes);
% loop on the elements
for i = 1:obj.mesh.number_elements
    % Element node numbers
    nodeA  =  obj.mesh.connect(i,2);
    nodeB  =  obj.mesh.connect(i,3);
    % index of nodes dof
    index_global_A = 3*(nodeA - 1) + [1:3];
    index_global_B = 3*(nodeB - 1) + [1:3];
    index = [index_global_A index_global_B];
    % element information (length and angle w.r.t horizontal)
    [L_element, theta_element] = obj.elementOrientation(i);
    % rotation matrix to go from local frame to global frame
    rot_matrix = obj.transformationMatrix(theta_element);
    % elementary matrices computation
    if strcmp(type,'mass')
        A_elem = obj.elementary_mass_matrix(L_element) ;
    elseif strcmp(type,'stiffness_at_origin')
        A_elem = obj.elementary_tangent_stiffness_matrix(L_element, zeros(6,1)); % elmentatry tangent stiffness at origin
    elseif strcmp(type,'stiffness_at_qs')
        qse = rot_matrix'*qs(index); % back to local coordinates
        A_elem = obj.elementary_tangent_stiffness_matrix(L_element, qse); % elmentatry tangent stiffness at qse
    else
        error('unknown matrix type')
    end
    % assemble the element matrix inside the full matrix
    A_full(index,index) = A_full(index,index) + rot_matrix*A_elem*rot_matrix'; % rotate back to global frame
end

end
