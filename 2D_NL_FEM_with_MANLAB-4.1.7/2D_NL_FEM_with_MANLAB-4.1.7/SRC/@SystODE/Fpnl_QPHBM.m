function [Fpnl_tot] = Fpnl_QPHBM(sys,U,p)
% Compute  the r.h.s Fpnl at order p from U (Taylor class)

H    = sys.H;
nz_tot   = sys.nz_tot;

coef_per_var = (H+1)*4*H+1;

Fpnl_fourier=zeros(coef_per_var,nz_tot);
aux_var = zeros(5,1);
cond_ini = zeros(1,nz_tot);
for r=1:p-1
    Ur  =get(U,'coefk',r) ;
    [Zr_tot,omega1r,omega2r,lambdar,omega1sqr,omega2sqr,omega1omega2r,lambdaomega1r,lambdaomega2r] = sys.get_Ztot(Ur);
    Upmr=get(U,'coefk',p-r);
    [Zpmr_tot,omega1pmr,omega2pmr] = sys.get_Ztot(Upmr);
    
    Zr_t0 = sum(Zr_tot(1:(coef_per_var-1)/2+1,:),1);
    Zpmr_t0 = sum(Zpmr_tot(1:(coef_per_var-1)/2+1,:),1);
    
    %boucle sur les listes
    for i=1:numel(sys.ic2)
        Fpnl_fourier(1,sys.ic2(i))= Fpnl_fourier(1,sys.ic2(i)) + sys.vc2(i)*lambdar*lambdapmr ;
    end
    
    for i=1:numel(sys.iforce2)
        Fpnl_fourier( (sys.hforce2(i)-1)*(2*H+1)+1,sys.iforce2(i))= Fpnl_fourier( (sys.hforce2(i)-1)*(2*H+1)+1,sys.iforce2(i)) - sys.vforce2(i)*lambdar*lambdapmr ;
    end
    
    for i=1:numel(sys.il1)
        Fpnl_fourier(:,sys.il1(i))= Fpnl_fourier(:,sys.il1(i)) + (sys.vl1(i)*lambdar)*Zpmr_tot(:,sys.jl1(i)) ;
    end
    
    for i=1:numel(sys.iq)
        Fpnl_fourier(:,sys.iq(i))= Fpnl_fourier(:,sys.iq(i)) + sys.vq(i)*( sys.Prod(Zr_tot(:,sys.jq(i)),Zpmr_tot(:,sys.kq(i))) ) ;
    end
    
    for i=1:numel(sys.id)
        Fpnl_fourier(:,sys.id(i))= Fpnl_fourier(:,sys.id(i)) + sys.vd(i)*sys.D(Zpmr_tot(:,sys.jd(i)),omega1r,omega2r) ;
    end
    
    for i=1:numel(sys.idd)
        Fpnl_fourier(:,sys.idd(i))= Fpnl_fourier(:,sys.idd(i)) + sys.vdd(i)*sys.DD(Zpmr_tot(:,sys.jdd(i)),omega1sqr,omega2sqr,omega1omega2r) ;
    end
    
    for i=1:numel(sys.id1)
        Fpnl_fourier(:,sys.id1(i))= Fpnl_fourier(:,sys.id1(i)) + sys.vd1(i)*sys.D(Zpmr_tot(:,sys.jd1(i)),lambdaomega1r,lambdaomega2r) ;
    end
    
    % definition of omega2 and lambdaomega :
    aux_var = aux_var - [omega1r*omega1pmr; omega2r*omega2pmr; omega1r*omega2pmr; lambdar*omega1pmr; lambdar*omega2pmr];
    
    if ~isempty(sys.idl)
        for i=1:numel(sys.idl)
            Fpnl_fourier(:,sys.idl(i))= Fpnl_fourier(:,sys.idl(i)) + sys.vdl(i)*sys.D(Zpmr_tot(:,sys.jdl(i)),omega1r,omega2r) ;
        end
        
        % Assumption : the variable after kdq(i) is the derivative of the variable kdq(i)
        for i=1:numel(sys.idq)
            Fpnl_fourier(:,sys.idq(i))= Fpnl_fourier(:,sys.idq(i)) + sys.vdq(i)*sys.Prod(Zr_tot(:,sys.jdq(i)),Zpmr_tot(:,sys.kdq(i)+1)) ;
        end
        
        cond_ini(sys.idq)= cond_ini(sys.idq) + (p-r)/p*sys.vdq'.*Zr_t0(sys.jdq).*Zpmr_t0(sys.kdq) ;
    end
end

% transform the matrix Fpnl to a column vector and add 0 for the phase
% equation, omega2,.. eqautions
Fpnl_fourier(1,sys.idq) = cond_ini(sys.idq);
Fpnl_fourier = reshape(Fpnl_fourier,nz_tot*coef_per_var,1);

Fpnl = [Fpnl_fourier(1:sys.neq-2);0;0];
Fpnl_aux = [Fpnl_fourier(sys.neq-1:end);aux_var];
Fpnl_tot = [Fpnl;Fpnl_aux];

