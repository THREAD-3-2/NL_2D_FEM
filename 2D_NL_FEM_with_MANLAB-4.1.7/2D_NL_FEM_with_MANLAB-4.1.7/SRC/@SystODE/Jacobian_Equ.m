function Jacobian = Jacobian_Equ(sys,Utot)
% Compute the Jacobian of R with respect to U

u = Utot(1:sys.nz);
Uaux = Utot(sys.nz+2:end);
Uvar = [u;Uaux];
lambda = Utot(sys.nz+1);

dRtotdZtot = sparse(sys.il0,sys.jl0,sys.vl0,sys.nz_tot,sys.nz_tot) ...
    + sparse(sys.il1,sys.jl1,sys.vl1*lambda,sys.nz_tot,sys.nz_tot) ...
    + sparse(sys.iq,sys.kq,2*sys.vq.*Uvar(sys.jq),sys.nz_tot,sys.nz_tot) ...
    + sparse(sys.idl,sys.jdl,sys.vdl,sys.nz_tot,sys.nz_tot) ...
    + sparse(sys.idq,sys.kdq,sys.vdq.*Uvar(sys.jdq),sys.nz_tot,sys.nz_tot) ;

dRtotdlambda = sparse(sys.il1,ones(1,size(sys.il1,1)),sys.vl1.*Uvar(sys.jl1),sys.nz_tot,1) ...
    + sparse(sys.ic1,ones(1,size(sys.ic1,1)),sys.vc1,sys.nz_tot,1) ...
    + sparse(sys.ic2,ones(1,size(sys.ic2,1)),(2*lambda)*sys.vc2,sys.nz_tot,1);

    Jacobian.dRdU = [dRtotdZtot(1:sys.nz,1:sys.nz) , dRtotdlambda(1:sys.nz)];
    Jacobian.dRdUaux = dRtotdZtot(1:sys.nz,sys.nz+1:end);
    Jacobian.dRauxdU = [dRtotdZtot(sys.nz+1:end,1:sys.nz) , dRtotdlambda(sys.nz+1:end)];
    Jacobian.dRauxdUaux = dRtotdZtot(sys.nz+1:end,sys.nz+1:end);
    Jacobian.dRtotdUtot = [Jacobian.dRdU , Jacobian.dRdUaux ; Jacobian.dRauxdU , Jacobian.dRauxdUaux];
    
end
