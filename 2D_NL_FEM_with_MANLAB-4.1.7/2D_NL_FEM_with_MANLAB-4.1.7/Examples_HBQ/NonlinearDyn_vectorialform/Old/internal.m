function [f_int_elem] = internal(params, c, s, eps, gam, k_bar, L_element)
f_int_elem = zeros(6,1);

f_int_elem(1) = -params.modulus*params.area*eps*c + ...
    params.shear*params.area*gam*s;
f_int_elem(2) = -params.modulus*params.area*eps*s + ...
   -params.shear*params.area*gam*c;
f_int_elem(3) =  params.modulus*params.area*eps*(L_element/2)*gam + ...
   -params.shear*params.area*(1 + eps)*(L_element/2)*gam + ...
   -params.modulus*params.I*k_bar;
f_int_elem(4) =  params.modulus*params.area*eps*c + ...
   -params.shear*params.area*gam*s;
f_int_elem(5) =  params.modulus*params.area*eps*s + ...
    params.shear*params.area*gam*c;
f_int_elem(6) =  params.modulus*params.area*eps*(L_element/2)*gam + ...
   -params.shear*params.area*(1 + eps)*(L_element/2)*gam + ...
    params.modulus*params.I*k_bar;
end