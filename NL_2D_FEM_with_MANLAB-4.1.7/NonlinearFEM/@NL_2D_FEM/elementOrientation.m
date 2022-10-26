function [L_element,theta_element] = elementOrientation(obj, i)
% compute the length and orientation angle of a given element number "i"
% gather the node numbers of the element
nodeA = obj.mesh.connect(i,2);
nodeB = obj.mesh.connect(i,3);
% define the element "geometric vector"
delta_x = obj.mesh.nodes(nodeB,2) - obj.mesh.nodes(nodeA,2);
delta_y = obj.mesh.nodes(nodeB,3) - obj.mesh.nodes(nodeA,3);
% delta_z = obj.mesh.nodes(nodeB,4) - obj.mesh.nodes(nodeA,4);
% length and orientation angle
% L_element = sqrt(delta_x^2 + delta_y^2 + delta_z^2);
L_element = sqrt(delta_x^2 + delta_y^2);
theta_element = atan2(delta_y,delta_x);
end

