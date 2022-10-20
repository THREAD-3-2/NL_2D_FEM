function v = atan(u)
% Taylor/exp est la formule de reu.coefrrence pour l'exponetielle
global Ck

g=1+u.*u;
atan0=atan(u.value);

atanu=u.coef*0; 
sto=1./g.value;
for k=1:Ck
    atanu(:,:,k)=0;
    for j=1:k-1
        atanu(:,:,k)=atanu(:,:,k)+j*times(atanu(:,:,j),g.coef(:,:,k-j));         
    end
    atanu(:,:,k)=sto.*(u.coef(:,:,k)-atanu(:,:,k)/k);
end
v = Taylor(u.order,atan0,atanu);