function [Fpnl_tot] = Fpnl_Equ(sys,U,p)
% Compute  the r.h.s Fpnl at order p from U (Taylor class)

Us = reshape(get(U,'coef'),sys.ninc,get(U,'order'));
Us = Us(:,1:p-1); % orders between 1 and p-1 are used only.

u = Us(1:sys.nz,:);
Uaux = Us(sys.nz+2:end,:);
Uvar = [u;Uaux];
lambda = Us(sys.nz+1,:);

% convolutions
valc2=lambda*lambda(p-1:-1:1)';
vall1=sum(lambda.*Uvar(sys.jl1,p-1:-1:1),2);

valq=sum(Uvar(sys.jq,:).*Uvar(sys.kq,p-1:-1:1),2);
valdq=sum(Uvar(sys.jdq,:).*Uvar(sys.kdq,p-1:-1:1).*(p-1:-1:1)/p,2);
                        
Fpnl_tot = sparse(sys.ic2,ones(1,size(sys.ic2,1)),sys.vc2*valc2,sys.neq_tot,1) ...
    + sparse(sys.il1,ones(1,size(sys.il1,1)),sys.vl1.*vall1,sys.neq_tot,1) ...
    + sparse(sys.iq,ones(1,size(sys.iq,1)),sys.vq.*valq,sys.neq_tot,1) ...
    + sparse(sys.idq,ones(1,size(sys.idq,1)),sys.vdq.*valdq,sys.neq_tot,1);

end