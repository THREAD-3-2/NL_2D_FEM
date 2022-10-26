function [val] = get(taylor,propName,index)
% GET Get TAYLOR properties from the specified object 
% and return the value
% Options are: order, value, coef, (coef1, coefk)
switch propName
    case 'order'
        val = taylor.order;
    case 'value'
        val = taylor.value;
    case 'coef'
        val = taylor.coef;
    case 'coef1'
        val = taylor.coef(:,:,1);
    case 'coefk'
        val = taylor.coef(:,:,index);
end