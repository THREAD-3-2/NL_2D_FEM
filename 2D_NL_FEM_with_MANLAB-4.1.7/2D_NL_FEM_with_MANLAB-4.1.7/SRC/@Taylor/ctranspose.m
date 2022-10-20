function v = ctranspose(u)
% @Taylor\ctranspose: {\tt u'} is the matrix transpose of {\tt u}. 
% For complex matrices, this is the complex conjugate transpose.
global Ck

cu = u.coef;
[m,n,p]=size(cu);
cv=zeros(n,m,p);
for k=1:Ck
    cv(:,:,k)=ctranspose(cu(:,:,k));
end
v = Taylor(u.order,ctranspose(u.value),cv);

