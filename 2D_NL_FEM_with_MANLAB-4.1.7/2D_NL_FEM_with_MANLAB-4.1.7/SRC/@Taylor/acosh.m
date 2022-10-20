function v = acosh(u)
% Taylor/exp est la formule de reu.coefrrence pour l'exponetielle
global Ck

g=sqrt(u*u-1);
acosh0=acosh(u.value);

acosh=u.coef*0; 
sto=1./g.value;
for k=1:Ck
    acosh(:,:,k)=0;
    for j=1:k-1
        acosh(:,:,k)=acosh(:,:,k)+j*times(acosh(:,:,j),g.coef(:,:,k-j));         
    end
    acosh(:,:,k)=sto*(u.coef(:,:,k)-acosh(:,:,k)/k);
end
v = Taylor(u.order,acosh0,acosh);