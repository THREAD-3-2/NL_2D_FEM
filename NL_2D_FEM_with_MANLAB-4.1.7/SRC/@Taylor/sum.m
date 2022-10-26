function r = sum(p)
% Taylor/sum: higher order recurrence formula sum(p)
global Ck

co = p.order;
cp = p.coef;
cp0 = p.value;

value=sum(p.value);
coef=zeros(1,1,co);
for k=1:Ck
    coef(1,1,k)=sum(p.coef(:,:,k));
end

r = Taylor(co,value,coef);
