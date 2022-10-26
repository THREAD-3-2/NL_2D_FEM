function [Fpnl_tot] = Fpnl(sys,U,p)
% Compute  the r.h.s Fpnl at order p from U (Taylor class)

% version 20/3/2017 
%
% Taylor coefficients of U until order p-1 are stored in Us
U1=get(U,'coef1');
Us=zeros(size(U1,1),p-1);
Us(:,1)=U1;
for r=2:p-1
 Us(:,r)= get(U,'coefk',r);
end

% Us = reshape(get(U,'coef'),obj.ninc,obj.order);
% Us = Us(:,1:p-1); % orders between 1 and p-1 are used only.

% convolution 
valQ=sum(Us(sys.jQ,:).*Us(sys.kQ,p-1:-1:1),2);
valdQ=sum(Us(sys.jdQ,:).*Us(sys.kdQ,p-1:-1:1).*(p-1:-1:1)/p,2);
                        
Fpnl_tot = sparse(sys.iQ,ones(1,size(sys.iQ,1)),sys.vQ.*valQ,sys.neq_tot,1) ...
    +  sparse(sys.idQ,ones(1,size(sys.idQ,1)),sys.vdQ.*valdQ,sys.neq_tot,1);

end