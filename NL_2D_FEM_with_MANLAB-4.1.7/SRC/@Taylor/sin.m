function v = sin(u)
% Taylor/exp est la formule de recurrence pour l'exponetielle
global Ck

co = u.order;
cu = u.coef; 
for k=1:Ck
    tcu(:,:,k)=k*cu(:,:,k);
end
cos0=cos(u.value);
sin0=sin(u.value);

c=u.coef*0; tc=c;
s=c; ts=s;
for k=1:Ck
    tc(:,:,k)=-times(tcu(:,:,k),sin0);
    ts(:,:,k)=times(tcu(:,:,k),cos0);
    for j=1:k-1
        tc(:,:,k)=tc(:,:,k)-times(tcu(:,:,j),s(:,:,k-j));         
        ts(:,:,k)=ts(:,:,k)+times(tcu(:,:,j),c(:,:,k-j)); 
    end
    c(:,:,k)=tc(:,:,k)/k;
    s(:,:,k)=ts(:,:,k)/k;
end
v = Taylor(co,sin0,s);