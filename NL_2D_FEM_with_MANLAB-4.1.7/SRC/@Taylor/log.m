function v = log(u)
% Taylor/log est la formule de recurrence pour le logarithme
global Ck

co = u.order;
cu = u.coef;

value=log(u.value);
du0=1./u.value;
cv=cu*0;
tcv=cv;
for k=1:Ck
    for j=1:k-1
        tcv(:,:,k) = tcv(:,:,k) + times(cu(:,:,k-j),tcv(:,:,j));  
    end
    tcv(:,:,k) = times(du0,(k*cu(:,:,k) - tcv(:,:,k)));
    cv(:,:,k) = tcv(:,:,k)/k;
end
v = Taylor(co,value,cv);



