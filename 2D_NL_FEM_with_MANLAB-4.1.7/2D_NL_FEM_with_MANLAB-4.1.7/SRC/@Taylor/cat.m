function v = cat(dim,u,w)
% Taylor/transpose: recurrence formula for higher order differentiation of
% p'
global Ck

if isa(u,'Taylor')
    co = u.order;
    u0=u.value;
    cu = u.coef;
    [mu,nu]=size(u0);
    if isa(w,'Taylor')
        w0=w.value;
        cw=w.coef;
        [mw,nw]=size(w0);
        if dim==1
            v0=zeros(mu+mw,nu);
            cv=zeros(mu+mw,nu,co);
        else
            v0=zeros(mu,nu+nw);
            cv=zeros(mu,nu+nw,co);
        end
        v0=cat(dim,u0,w0);
        for k=1:Ck
            cv(:,:,k)=cat(dim,cu(:,:,k),cw(:,:,k));
        end
    else
       [mw,nw]=size(w);
        if dim==1
            v0=zeros(mu+mw,nu);
            cv=zeros(mu+mw,nu,co);
        else
            v0=zeros(mu,nu+nw);
            cv=zeros(mu,nu+nw,co);
        end
        v0=cat(dim,u0,w);
        for k=1:Ck
            cv(:,:,k)=cat(dim,cu(:,:,k),zeros(mw,nw));
        end
    end
else
    co=w.order;
    w0=w.value;
    cw=w.coef;
    [mu,nu]=size(u);
    [mw,nw]=size(w0); 
    if dim==1
        v0=zeros(mu+mw,nu);
        cv=zeros(mu+mw,nu,Ck);
    else
        v0=zeros(mu,nu+nw);
        cv=zeros(mu,nu+nw,Ck);
    end
    v0=cat(dim,u,w0);
    for k=1:Ck
        cv(:,:,k)=cat(dim,zeros(mu,nu),cw(:,:,k));
    end 
end
v = Taylor(co,v0,cv);
