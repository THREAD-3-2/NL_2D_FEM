function r = mpower(p,a)
% Taylor/mpower: higher order recurrence formula for p^a(Leibniz)
global Ck

if (a == 2)
    r=p;
    cp = p.coef;
    p0 = p.value;
    r.value = p0*p0;
    cr=cp;
    for k=1:Ck
        cr(:,:,k) = 2*p0(:,:)*cp(:,:,k);
    end
    for k=1:Ck
        for j = 1:k-1
            cr(:,:,k) = cr(:,:,k)+cp(:,:,j)*cp(:,:,k-j);  
        end
    end
    r.coef=cr;
else
    if isa(a,int8)
        r=p*p^(a-1);
    else
        fprintf('Taylor/mpower not yet implemented for exponent~=2')
        stop
    end
end
