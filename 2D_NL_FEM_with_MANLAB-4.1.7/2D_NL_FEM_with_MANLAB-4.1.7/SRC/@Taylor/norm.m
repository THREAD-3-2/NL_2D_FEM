function r = norm(p)
% Taylor/norm: recurrence formula for higher order differentiation of
% norm(p)
global Ck

co = p.order;
cp = p.coef;
cu=zeros(1,1,co);
cv=cu;
%differentiation de dot
u0 = dot(p.value,p.value);
for k=1:Ck
    cu(1,1,k)= 2*dot(p.value,cp(:,:,k));
    for l=1:k-1
        cu(1,1,k)=cu(1,1,k) + dot(cp(:,:,l),cp(:,:,k-l));
    end
end
%differentation de sqrt
v0=sqrt(u0);
if v0~=0
    r0=0.5/v0;
    for k=1:Ck
        cv(1,1,k)=cu(1,1,k);
        for j=1:k-1
            cv(1,1,k)=cv(1,1,k)-cv(1,1,j)*cv(1,1,k-j);
        end
        cv(1,1,k)=cv(1,1,k)*r0;
    end
end

%fprintf('verifier norm dans @Taylor \n');
r = Taylor(co,v0,cv);

