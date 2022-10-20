function r=power(x,a)
%MPOWER Elementwise power (.^) operator for Taylor objects.
%
%LAMPOH K. 2013

global Ck

if isa(x,'Taylor')
    if isa(a,'Taylor')
        error('This mpower is not already implemented, the exponent must be double')
    else
        r.value=x.value.^a;

        r.coef=zeros([size(r.value),x.order]);
        r.order=x.order;
        r.coef(:,:,1)=a*r.value.*x.coef(:,:,1)./x.value;

        for kk=2:Ck
            for jj=1:kk-1
                r.coef(:,:,kk)=r.coef(:,:,kk)+(a*(kk-jj)-jj)...
                    *r.coef(:,:,jj).*x.coef(:,:,kk-jj);
            end
            r.coef(:,:,kk)=(r.coef(:,:,kk)+a*kk...
                *r.value.*x.coef(:,:,kk))./(kk*x.value);
        end
    end
end

r=Taylor(x.order,r.value,r.coef);
