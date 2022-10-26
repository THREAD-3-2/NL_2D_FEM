function [val] = evalseries(t,a,k)
% Evaluate the taylor series stored it t at point a, up to order k

val = t.value; ak = 1;
for l=1:k
    ak = ak*a;
	val = val + t.coef(:,:,l)*ak;
end

