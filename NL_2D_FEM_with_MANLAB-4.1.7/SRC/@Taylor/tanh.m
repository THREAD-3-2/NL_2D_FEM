function vh = tanh(u)
% Taylor/exp est la formule de reu.coefrrence pour l'exponetielle
global Ck

gh=cosh(u).*cosh(u);
tanh0=tanh(u.value);

thu=u.coef*0; 
stoh=1./gh.value;
for k=1:Ck
    thu(:,:,k)=0;
    for j=1:k-1
        thu(:,:,k)=thu(:,:,k)+j*times(thu(:,:,j),gh.coef(:,:,k-j));         
    end
    thu(:,:,k)=stoh.*(u.coef(:,:,k)-thu(:,:,k)/k);
end
vh = Taylor(u.order,tanh0,thu);