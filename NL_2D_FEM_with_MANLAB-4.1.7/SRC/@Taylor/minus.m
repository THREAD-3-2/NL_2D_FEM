function r = minus(p,q)
% Taylor/minus: recurrence formula for higher order differentiation of p-q
% Matlab definition:
% Subtraction. A-B subtracts B from A. A and B must have the same size,
% unless one is a scalar. A scalar can be subtracted from a matrix of any size
global Ck

if isa(p,'Taylor')
    r=p;
    p0=p.value;
    if isa(q,'Taylor')
        r.value=p0-q.value;
        r.coef(:,:,1:Ck)=p.coef(:,:,1:Ck)-q.coef(:,:,1:Ck);
        
    else 
        [mp,np]=size(p0);
        [mq,nq]=size(q);
        
        if ([mp,np] == [mq,nq])
            r.value=p0-q;
        else
            if ([mq,nq] == [1,1])
                r.value=p0-q*ones(mp,np);
            else
                fprintf('\n @Taylor:moins dimension does not match\n')
            end
        end
    end
else
    r=-q;
    q0=r.value;
    [mp,np]=size(p);
    [mq,nq]=size(q0);
    if ([mp,np] == [mq,nq])
        r.value=p+q0;
    else
        if ([mp,np] == [1,1])
            r.value=p*ones(mq,nq)+q0;
        else
            fprintf('\n @Taylor:moins dimension does not match\n')
        end
    end
end