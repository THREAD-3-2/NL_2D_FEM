function v = acos(u)
% Taylor/exp est la formule de reu.coefrrence pour l'exponetielle
global Ck

g=sqrt(1-u.*u); 
acos0=acos(u.value);

acosu=u.coef*0; 
sto=-1./g.value;
for k=1:Ck
    acosu(:,:,k)=0;
    for j=1:k-1
        acosu(:,:,k)=acosu(:,:,k)+j*times(acosu(:,:,j),g.coef(:,:,k-j));         
    end
    acosu(:,:,k)=sto.*(u.coef(:,:,k)+acosu(:,:,k)/k);
end
v = Taylor(u.order,acos0,acosu);