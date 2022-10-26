function r = times(p,q)
% Taylor/mtimes: higher order recurrence formula for p*q (Leibniz)
global Ck

if isa(p,'Taylor')
    co = p.order;
    cp = p.coef;
    cp0 = p.value;
    %[m,n,l]=size(p.coef);
    %cr=zeros(m,n,l);
    cr=cp*0;
    if isa(q,'Taylor')
        cq0 = q.value;
        cq = q.coef; 
        value = times(cp0,cq0);
        for k=1:Ck
            cr(:,:,k) = times(cp0,cq(:,:,k)) + times(cp(:,:,k),cq0);
            for j = 1:k-1
                cr(:,:,k) = cr(:,:,k)+times(cq(:,:,j),cp(:,:,k-j));  
            end
        end
    else 
        value=times(cp0,q);
        for k=1:Ck
            cr(:,:,k)=times(cp(:,:,k),q);
        end
    end
else
    co = q.order;
    cq = q.coef;
    value=times(p,q.value);
    cr=cq*0;
    for k=1:Ck
        cr(:,:,k)=times(p,cq(:,:,k));
    end
end
r = Taylor(co,value,cr);
%fprintf('mtimes: a verifier \n')
