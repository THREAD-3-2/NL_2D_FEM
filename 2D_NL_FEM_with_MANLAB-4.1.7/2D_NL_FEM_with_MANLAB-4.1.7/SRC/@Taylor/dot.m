function r = dot(p,q)
% Taylor/dot: higher order recurrence formula for dot(p,q)
global Ck

if isa(p,'Taylor')
    co = p.order;
    cp = p.coef; 
    cr=zeros(1,1,co);
    if isa(q,'Taylor') 
        cq = q.coef; 
        value = dot(p.value,q.value);
        for k=1:Ck
            cr(1,1,k)= dot(p.value,cq(:,:,k))+dot(cp(:,:,k),q.value);
            for l=1:k-1
                cr(1,1,k)=cr(1,1,k) + dot(cp(:,:,l),cq(:,:,k-l));
            end
        end
    else
        cr=zeros(1,1,co);
        value=dot(p.value(:),q);
        for k=1:Ck
            cr(1,1,k)=dot(cp(:,:,k),q);
        end
    end
else
    co = q.order;
    cq = q.coef; 
    cr=zeros(1,1,co);
    value=dot(p,q.value(:));
    for k=1:Ck
        cr(1,1,k)=dot(cq(:,:,k),p);
    end
end
%fprintf('     dot\n')
r = Taylor(co,value,cr);
