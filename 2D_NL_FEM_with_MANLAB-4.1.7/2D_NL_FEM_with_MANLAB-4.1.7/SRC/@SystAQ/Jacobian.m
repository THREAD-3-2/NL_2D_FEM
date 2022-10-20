function [Jacobian] = Jacobian(sys,Utot)
% Compute the Jacobian of R with respect to Utot
%   dRtotdUtot =  L(.)   + 2*Q(U, . )  + dL(.) + dQ(U,.)
% It gives the block decomposition of the jacobian matrix in the field of
% the structure "Jacobian".

dRtotdUtot = sparse(sys.iL,sys.jL,sys.vL,sys.neq_tot,sys.ninc) ...
    + sparse(sys.iQ,sys.kQ,2*sys.vQ.*Utot(sys.jQ),sys.neq_tot,sys.ninc) ...
    + sparse(sys.idL,sys.jdL,sys.vdL,sys.neq_tot,sys.ninc) ...
    + sparse(sys.idQ,sys.kdQ,sys.vdQ.*Utot(sys.jdQ),sys.neq_tot,sys.ninc) ;

    Jacobian.dRtotdUtot = dRtotdUtot;
    Jacobian.dRdU = dRtotdUtot(1:sys.neq,1:sys.neq+1);
    Jacobian.dRdUaux = dRtotdUtot(1:sys.neq,sys.neq+2:end);
    Jacobian.dRauxdU = dRtotdUtot(sys.neq+1:end,1:sys.neq+1);
    Jacobian.dRauxdUaux = dRtotdUtot(sys.neq+1:end,sys.neq+2:end);

end


