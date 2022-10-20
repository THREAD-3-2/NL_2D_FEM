function [val] = evalderiv(t,a,k)
% Evaluate the derivative of taylor series,  at point a, up to order k

val = t.value*0 ; ak = 1;
for l=1:k 
       val = val + t.coef(:,:,l)*ak;
        ak = ak*a*(l+1)/l;    
end

