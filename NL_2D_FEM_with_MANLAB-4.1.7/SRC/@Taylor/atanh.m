function v = atanh(u)
% Taylor/exp est la formule de reu.coefrrence pour l'exponetielle
global Ck

g=1-u.*u;
atanh0=atanh(u.value);

atanh=u.coef*0; 
sto=1./g.value;
for k=1:Ck
    atanh(:,:,k)=0;
    for j=1:k-1
        atanh(:,:,k)=atanh(:,:,k)+j*times(atanh(:,:,j),g.coef(:,:,k-j));         
    end
    atanh(:,:,k)=sto.*(u.coef(:,:,k)-atanh(:,:,k)/k);
end
v = Taylor(u.order,atanh0,atanh);