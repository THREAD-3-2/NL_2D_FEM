function [Res, Jac] = man_residual_fixed_frequency(U0, omega0, sys)
% residual for fixed frequancy MAN
Res = [sys.R(sys, U0); U0((2*sys.H + 1)*sys.nz+1)-omega0];

if nargout>1
J = sys.Jacobian(U0).dRtotdUtot;
en = zeros(1, size(J,2));
en( (2*sys.H + 1)*sys.nz+1 ) = 1;
Jac = [ J;
        en];
end
end