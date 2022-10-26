function [T] = transformationMatrix(obj, theta_element)
% transformationMatrix Calculates the rotation transformation matrix 
% [T] for a given element orientation angle

R = [ cos(theta_element) -sin(theta_element)     0
      sin(theta_element)  cos(theta_element)     0
               0                     0           1 ];
N = zeros(3,3);
           
T = [ R N;
      N R ];
end

