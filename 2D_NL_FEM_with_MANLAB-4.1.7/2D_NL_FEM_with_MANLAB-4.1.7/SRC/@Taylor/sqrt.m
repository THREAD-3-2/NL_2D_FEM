function v = sqrt(u)
% Taylor/mrdivide: higher order recurrence formula for p/q

global Ck

cu(:,:,:) = u.coef(:,:,:);
cv(:,:,:) = 0*cu(:,:,:);
v0=sqrt(u.value);
dw0 = 1./(2*v0);    %initialisation

% cv(:,:,1:Ck)=cu(:,:,1:Ck); deja fait avec v=u
for k=1:Ck
    cv(:,:,k)=cu(:,:,k);
    for j=1:k-1
        cv(:,:,k)=cv(:,:,k)-cv(:,:,j)*cv(:,:,k-j);
    end
    cv(:,:,k) = dw0*cv(:,:,k);
end
v=Taylor(u.order,v0,cv);
