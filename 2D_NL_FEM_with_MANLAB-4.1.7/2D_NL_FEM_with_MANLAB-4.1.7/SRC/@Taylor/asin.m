function v = asin(u)
% Taylor/exp est la formule de reu.coefrrence pour l'exponetielle
global Ck

g=sqrt(1-u*u);
asin0=asin(u.value);

asinu=u.coef*0; 
sto=1./g.value;
for k=1:Ck
    asinu(:,:,k)=0;
    for j=1:k-1
        asinu(:,:,k)=asinu(:,:,k)+j*times(asinu(:,:,j),g.coef(:,:,k-j));         
    end
    asinu(:,:,k)=sto*(u.coef(:,:,k)-asinu(:,:,k)/k);
end
v = Taylor(u.order,asin0,asinu);