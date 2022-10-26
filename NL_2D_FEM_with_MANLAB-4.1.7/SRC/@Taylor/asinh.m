function v = asinh(u)
% Taylor/exp est la formule de reu.coefrrence pour l'exponetielle
global Ck

g=sqrt(1+u*u);
asinh0=asinh(u.value);

asinh=u.coef*0; 
sto=1./g.value;
for k=1:Ck
    asinh(:,:,k)=0;
    for j=1:k-1
        asinh(:,:,k)=asinh(:,:,k)+j*times(asinh(:,:,j),g.coef(:,:,k-j));         
    end
    asinh(:,:,k)=sto*(u.coef(:,:,k)-asinh(:,:,k)/k);
end
v = Taylor(u.order,asinh0,asinh);