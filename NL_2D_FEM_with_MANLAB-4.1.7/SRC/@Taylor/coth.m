function vh = coth(u)
% Taylor/exp est la formule de reu.coefrrence pour l'exponetielle
global Ck

gh=-sinh(u).*sinh(u);
coth0=coth(u.value);

thu=u.coef*0; 
stoh=1./gh.value;
for k=1:Ck
    thu(:,:,k)=0;
    for j=1:k-1
        thu(:,:,k)=thu(:,:,k)+j*times(thu(:,:,j),gh.coef(:,:,k-j));         
    end
    thu(:,:,k)=-stoh.*(u.coef(:,:,k)+thu(:,:,k)/k);
end
vh = Taylor(u.order,coth0,thu);