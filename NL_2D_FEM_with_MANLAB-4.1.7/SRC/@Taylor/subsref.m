function a = subsref(t,s)

switch s(1).type
    case '.'
        fieldname = s(1).subs;
        switch fieldname
            case 'value'
                a = cr_subsubsref(t.value,s);
                return;
            case 'coef'
                a = cr_subsubsref(t.coef,s);
                return;
            case 'order'
                a = cr_subsubsref(t.order,s);
                return;
            otherwise
                error(['No such field: ', fieldname]);
        end
end

a.value = t.value(s.subs{:});
a.coef=zeros([size(a.value),t.order]);

%for
    kk=1:t.order;
    %acoefkk=t.coef(:,:,kk);
    a.coef(:,:,kk)= t.coef(s.subs{:},kk);
%end
a=Taylor(t.order,a.value,a.coef);



