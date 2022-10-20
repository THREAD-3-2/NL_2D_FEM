% Computes the tangent vector for the point
% corresponding to the path parameter value  'a'
% using the power series

function [Utb]=tangenteserie(U,a)

Utb = zeros(size(U.value));
Norder=get(U,'order');
for i = 1:Norder-1
  Utb = Utb + i * a^(i-1) * U.coef(:,:,i+1);
end
