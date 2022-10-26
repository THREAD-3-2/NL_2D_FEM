
function r = subsasgn(r,s,val)
   % Implement a special subscripted assignment
  
switch s(1).type
case '()'
    ind = s.subs{:};
    if isa(val,'Taylor')         
        if isa(r, 'Taylor')           
            co=r.order;
            value=r.value;
            coef=r.coef;     
            value(ind,:,:)=val.value;    
            coef(ind,:,:)=val.coef;
            r=Taylor(co,value,coef);
        else
            co=val.order;
            value(ind,:,:)=val.value;
            coef(ind,:,:)=val.coef;         
            r=Taylor(co,value,coef);           
        end
    end
                  
    case '.'
        fprintf(' subsasgn . non traite dans Taylor')
        switch s(1).subs
            case 'coef'
                r = val.coef;
            case 'value'
                r=val.value;
            case 'order'
                r=order;
            otherwise
                error('erreur dans Taylor/subsasgn')
        end
end