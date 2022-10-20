function v = exp(u)
% Taylor/exp est la formule de recurrence pour l'exponetielle
global Ck

co = u.order;
cu = u.coef; 
for k=1:Ck
    tcu(:,:,k)=k*cu(:,:,k);
end
value=exp(u.value);
cv=cu*0; tcv=cv;
for k=1:Ck
    tcv(:,:,k) = times(value,tcu(:,:,k));
    for j=1:k-1
        tcv(:,:,k) = tcv(:,:,k) + times(cv(:,:,k-j),tcu(:,:,j));  
    end
    cv(:,:,k) = tcv(:,:,k)/k;
end
v = Taylor(co,value,cv);
