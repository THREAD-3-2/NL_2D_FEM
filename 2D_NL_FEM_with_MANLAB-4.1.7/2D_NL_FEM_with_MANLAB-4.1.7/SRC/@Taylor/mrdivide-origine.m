function v = mrdivide(u,w)
% Taylor/mrdivide: higher order recurrence formula for p/q

global Ck
if isa(w,'Taylor')
    v=w;
    cw(:,:,:) = w.coef(:,:,:);
    cv=cw*0;
    dw0 = 1./w.value;    %initialisation
    if isa(u,'Taylor') 
        v0(:,:) = times(u.value(:,:),dw0); v.value=v0;
        cu(:,:,:) = u.coef(:,:,:);
        for k=1:Ck
            cv(:,:,k)=times(v0(:,:),cw(:,:,k)); %j=0 de la formule
            for j=1:k-1
                cv(:,:,k)=cv(:,:,k)+times(cv(:,:,j),cw(:,:,k-j));
            end
            cv(:,:,k) = times(dw0,(u.coef(:,:,k)-cv(:,:,k)));
        end
        v.coef=cv;
    else
        v=w;
        cw(:,:,:) = w.coef(:,:,:);
        dw0 = 1./w.value;    %initialisation 
        v.value(:,:) = dw0;
        v0=dw0;
        for k=1:Ck
            cv(:,:,k)=times(v0(:,:),cw(:,:,k)); %j=0 de la formule
            for j=1:k-1
                cv(:,:,k)=cv(:,:,k)+times(cv(:,:,j),cw(:,:,k-j));
            end
            cv(:,:,k) = times(dw0,-cv(:,:,k));
        end
        v.coef=cv;
    end
else
    sto = 1./w;
    v=u;
    v.value = u.value*sto;
    v.coef=u.coef*sto;
end
