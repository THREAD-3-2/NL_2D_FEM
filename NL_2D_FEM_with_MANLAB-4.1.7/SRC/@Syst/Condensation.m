function [U,Uaux] = Condensation(sys,Jacobian,B,Ba,val)
%function [U,Uaux] = Condensation(sys,Jacobian,B,Ba,val)
%   Give Utot such that Jacobian*Utot = Btot and Jacobian a structure containing the parameters
%   of the condensation (LU decomposition of the matrices to
%   inverse,...).
%   Jacobian is a structure that contain either the sub-matrices
%   and the condensation informations.
%   If a value for lambda is specified, it is added at the end of the main
%   variables of Utot.

if (nargin == 3) || isempty(Jacobian.dRauxdUaux)
    % Here there are no auxiliary variables, we use standard formulaes to
    % resolve the linear system :
    if nargin == 4
        u(Jacobian.qK,1) = Jacobian.UK\(Jacobian.LK\B(Jacobian.pK));
        U=u;
    else
        u(Jacobian.qK,1) = Jacobian.UK\(Jacobian.LK\(B(Jacobian.pK) - val*Jacobian.missingcolumn(Jacobian.pK)));
        U = zeros(sys.neq+1,1);
        U(Jacobian.index_exchange) = val;
        U(Jacobian.list_val) = u;
    end
    
    Uaux=[];
    
else % Condensation of the auxiliary variables
    
    if nargin == 4
        % Right-hand-side / Condensation formulaes
        Rhs = B - Jacobian.dRdUaux(:,Jacobian.qK_aux)*(Jacobian.UK_aux\(Jacobian.LK_aux\Ba(Jacobian.pK_aux)));
        
        % Computation of the main variables, then the auxiliary variables :
        u(Jacobian.qK,1) = Jacobian.UK\(Jacobian.LK\Rhs(Jacobian.pK));
        Uaux(Jacobian.qK_aux,1) = Jacobian.UK_aux\(Jacobian.LK_aux\(Ba(Jacobian.pK_aux) - Jacobian.dRauxdU(Jacobian.pK_aux,:)*u));
        U = u;

    else
        % Right-hand-side / Condensation formulaes
        Rhs = B - Jacobian.dRdUaux(:,Jacobian.qK_aux)*(Jacobian.UK_aux\(Jacobian.LK_aux\Ba(Jacobian.pK_aux))) - val*Jacobian.missingcolumn;
        
        % Computation of the main variables, then the auxiliary variables :
        u(Jacobian.qK,1) = Jacobian.UK\(Jacobian.LK\Rhs(Jacobian.pK));
        U = zeros(sys.neq+1,1);
        U(Jacobian.index_exchange) = val;
        U(Jacobian.list_val) = u;
        Uaux(Jacobian.qK_aux,1) = Jacobian.UK_aux\(Jacobian.LK_aux\(Ba(Jacobian.pK_aux) - Jacobian.dRauxdU(Jacobian.pK_aux,:)*U));
    end
end


end

