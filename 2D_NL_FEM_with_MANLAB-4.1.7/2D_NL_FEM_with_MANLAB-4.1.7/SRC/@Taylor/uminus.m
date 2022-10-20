function r = uminus(p)
% Taylor/minus: recurrence formula for higher order differentiation of -p
global Ck

r = p;
r.value = -r.value;
r.coef(:,:,1:Ck)=-r.coef(:,:,1:Ck);
