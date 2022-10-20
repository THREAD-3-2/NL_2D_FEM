function r = mtimes(p,q)
% Taylor/mtimes: higher order recurrence formula for p*q (Leibniz)
global Ck

if isa(p,'Taylor')
    co = p.order;
    cp = p.coef;
    p0 = p.value;
    if isa(q,'Taylor')
        q0 = q.value;
        cq = q.coef; 
        value = p0*q0;
        [mp,nq]=size(value);
        cr=zeros(mp,nq,co);
        for k=1:Ck
            cr(:,:,k) = p0(:,:)*cq(:,:,k) + cp(:,:,k)*q0(:,:);
            for j = 1:k-1
                cr(:,:,k) = cr(:,:,k)+cp(:,:,k-j)*cq(:,:,j);  
            end
        end
        r = Taylor(co,value,cr);
    else %( p is Taylor, q is double)
        value=p0*q;
        [m,n]=size(value);
        r=Taylor(co,zeros(m,n));
        r.value=value;
        for k=1:Ck
            r.coef(:,:,k)=cp(:,:,k)*q;
        end
    end
else
    co = q.order;
    value=p*q.value;
    [m,n]=size(value);
    r=Taylor(co,zeros(m,n));
    r.value=value;
    for k=1:Ck
        r.coef(:,:,k)=p*q.coef(:,:,k);
    end
end