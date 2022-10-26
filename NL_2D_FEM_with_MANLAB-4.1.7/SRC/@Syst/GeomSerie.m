function [Uclean,Uscale, Alpha, status] = GeomSerie(sys, U)
% look for an emerging geometrical series in the last term of U (Taylor)
% If true, the geometrical serie is extracted
% outputs : Uclean (cleaned series), Uscale (scale factor),
%           Alpha (common ratio), status (flag)
Norder=sys.order;
Uclean=U;
Ucoef=get(U,'coef');
e_gp1 = 1e-5; %
e_gp2 = 1e-4; %

SumRatioUp=0;
alphap=zeros(3,1);

if (Norder >=10)
    for p=1:3
        alphap(p)=Ucoef(:,:,Norder-p)'*Ucoef(:,:,Norder)/(Ucoef(:,:,Norder)'*Ucoef(:,:,Norder));
        UpOrt=Ucoef(:,:,Norder-p)-alphap(p)*Ucoef(:,:,Norder);    % remainder
        ratioUp = norm(UpOrt) / norm(Ucoef(:,:,Norder-p));
        SumRatioUp = SumRatioUp  + ratioUp;
    end
    alphap=abs(alphap);
    EcarAlphap= ((alphap(2)^(1/2)/alphap(1)- 1))^2 +((alphap(3)^(1/3)/alphap(1)- 1))^2;
    
    % Dectection of geometric progression test
    if (EcarAlphap <= e_gp1 && SumRatioUp <= e_gp2)
        status='simplebif';
        Alpha=Ucoef(:,:,Norder)'*Ucoef(:,:,Norder-1)/(Ucoef(:,:,Norder)'*Ucoef(:,:,Norder));
        Uscale=Ucoef(:,:,Norder)*Alpha^Norder;
        
        disp([' Emerging geometrical serie detected, commun ratio: Alpha= ' num2str(Alpha)])   %
        % Extraction and suppession of the geometric serie from U.
        for k=1:Norder
            UCoefk = Ucoef(:,:,k)- Uscale/(Alpha^k);
            Uclean=set(Uclean,'coefk',UCoefk,k);
        end
    else
        status='nothing';
        Uscale=0;
        Alpha=0;
    end
    
else
    status='nothing';
    Uscale=0;
    Alpha=0;
end
end

