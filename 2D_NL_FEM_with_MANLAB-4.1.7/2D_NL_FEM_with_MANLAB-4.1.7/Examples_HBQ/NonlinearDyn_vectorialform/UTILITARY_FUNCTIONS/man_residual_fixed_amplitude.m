function [Res, Jac] = man_residual_fixed_amplitude(U0, idx, amp, sys)
% residual for fixad amplitude MAN equation
Res = [sys.R(sys, U0); U0(idx)-amp];
% Jacobian if needed
if nargout>1
    J = sys.Jacobian(U0).dRtotdUtot;
    en = zeros(1, size(J,2));
    en( idx ) = 1;
    Jac = [ J;
        en];
end
end