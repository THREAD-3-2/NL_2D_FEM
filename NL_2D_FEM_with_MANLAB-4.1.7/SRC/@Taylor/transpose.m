function v = transpose(u)
% @Taylor\transpose: {\tt u'} is the array transpose of {\tt u}. 
% For complex matrices, this does not involve conjugation.
global Ck

cu = u.coef;
[m,n,p]=size(cu);
cv=zeros(n,m,p);
for k=1:Ck
    cv(:,:,k)=transpose(cu(:,:,k));
end
v = Taylor(u.order,transpose(u.value),cv);

