function vh = cosh(u)
% Taylor/exp est la formule de recurrence pour l'exponetielle
global Ck

co = u.order;
cu = u.coef; 
for k=1:Ck
    tcu(:,:,k)=k*cu(:,:,k);
end
cosh0=cosh(u.value);
sinh0=sinh(u.value);

ch=u.coef*0; tch=ch;
sh=ch; tsh=sh;
for k=1:Ck
    tch(:,:,k)=times(tcu(:,:,k),sinh0);
    tsh(:,:,k)=times(tcu(:,:,k),cosh0);
    for j=1:k-1
        tch(:,:,k)=tch(:,:,k)+times(tcu(:,:,j),sh(:,:,k-j));         
        tsh(:,:,k)=tsh(:,:,k)+times(tcu(:,:,j),ch(:,:,k-j)); 
    end
    ch(:,:,k)=tch(:,:,k)/k;
    sh(:,:,k)=tsh(:,:,k)/k;
end
vh = Taylor(co,cosh0,ch);