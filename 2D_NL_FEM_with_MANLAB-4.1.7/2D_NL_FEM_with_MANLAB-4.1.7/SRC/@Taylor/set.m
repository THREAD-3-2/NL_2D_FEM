function taylor = set(taylor,prop,val,index)
% SET Set taylor properties and return the updated object
% Options are: order, value, coef, (coef1 for the first order information)
switch prop
case 'order'
    taylor.order = val;
case 'value'
    taylor.value = val;
case 'coef'
    taylor.coef = val;
case 'coef0'
    taylor.coef(:,:,:) = 0;
case 'coef1'
    taylor.coef(:,:,1) = val;
case 'coefk'
    taylor.coef(:,:,index) = val;
end