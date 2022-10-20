function [U,Uaux] = CondensationBif(sys,Jacobian,B,Ba)
%function [U,Uaux] = Condensation(sys,Jacobian,B,Ba,val)
%   Give Utot such that Jacobian*Utot = Btot and Jacobian a structure containing the parameters
%   of the condensation (LU decomposition of the matrices to
%   inverse,...).
%   Jacobian is a structure that contain either the sub-matrices
%   and the condensation informations.
%   If a value for lambda is specified, it is added at the end of the main
%   variables of Utot.


if  isempty(Jacobian.dRauxdUaux)
    % Here there are no auxiliary variables, we use standard formulaes to
    % resolve the linear system :
    U(Jacobian.qK,1) = Jacobian.UK\(Jacobian.LK\B(Jacobian.pK));    
    Uaux=[];
    
else % Condensation of the auxiliary variables
    
        % Right-hand-side / Condensation formulaes
        if ~istril(Jacobian.dRauxdUaux)
         Rhs = B - Jacobian.dRdUaux(:,Jacobian.qK_aux)*(Jacobian.UK_aux\(Jacobian.LK_aux\Ba(Jacobian.pK_aux)));
        else
         Rhs = B - Jacobian.dRdUaux*(Jacobian.dRauxdUaux\Ba) ;    
        end
   
        % Computation of the main variables, then the auxiliary variables :
        u = zeros(sys.neq,1);
        u(Jacobian.qK) = Jacobian.UK\(Jacobian.LK\Rhs(Jacobian.pK));
        Uaux = zeros(sys.neq_aux,1);
        Rhs_aux=Ba - Jacobian.dRauxdU*u;
        if ~istril(Jacobian.dRauxdUaux)
          Uaux(Jacobian.qK_aux) = Jacobian.UK_aux\(Jacobian.LK_aux\Rhs_aux(Jacobian.pK_aux)); 
        else
          Uaux = Jacobian.dRauxdUaux\Rhs_aux;   
        end
        U = u;
end


end

