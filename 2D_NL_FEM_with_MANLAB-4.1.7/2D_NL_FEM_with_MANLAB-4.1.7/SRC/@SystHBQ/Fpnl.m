function [Fpnl_tot] = Fpnl(sys,U,p)
% Compute  the r.h.s Fpnl at order p from U (Taylor class)

H    = sys.H;
nz_tot   = sys.nz_tot;

DHp1=2*H+1 ;

Fpnl_fourier=zeros(DHp1,nz_tot);
Fpnl_phase = 0;
aux_var = zeros(2,1);
cond_ini = zeros(1,nz_tot);

%%% As Q is symmetric, the sum is separated in two parts to optimize the
%%% time-computation.
for r=1:floor((p-1)/2)
    Ur  =get(U,'coefk',r) ;
    [Zr_tot,omegar,lambdar,omega2r,lambdaomegar] = sys.get_Ztot(Ur);
    Upmr=get(U,'coefk',p-r);
    [Zpmr_tot,omegapmr,lambdapmr,omega2pmr,lambdaomegapmr] = sys.get_Ztot(Upmr);
    
    % Vectorial computation of all the complex forms.
    Zr_compl = sys.real_to_compl(Zr_tot,'full');
    Zpmr_compl = sys.real_to_compl(Zpmr_tot,'full');
    if ~isempty(sys.id) || ~isempty(sys.idd) || ~isempty(sys.id1) ; Zd_pmr = sys.D(Zpmr_tot); Zd_r = sys.D(Zr_tot); end
    if ~isempty(sys.idd) ; Zdd_pmr = sys.D(Zd_pmr); Zdd_r = sys.D(Zd_r); end
    
    Zr_t0 = sum(Zr_tot(1:H+1,:),1);
    Zpmr_t0 = sum(Zpmr_tot(1:H+1,:),1);
    
    % Loops on the lists.
    for i=1:numel(sys.iforce2)
        Fpnl_fourier(sys.hforce2(i),sys.iforce2(i))= Fpnl_fourier(sys.hforce2(i),sys.iforce2(i)) - 2*sys.vforce2(i)*lambdar*lambdapmr ;
    end
    
    for i=1:numel(sys.ic2)
        Fpnl_fourier(1,sys.ic2(i))= Fpnl_fourier(1,sys.ic2(i)) + 2*sys.vc2(i)*lambdar*lambdapmr ;
    end
    
    for i=1:numel(sys.il1)
        Fpnl_fourier(:,sys.il1(i))= Fpnl_fourier(:,sys.il1(i)) + sys.vl1(i)*(lambdar*Zpmr_tot(:,sys.jl1(i)) + lambdapmr*Zr_tot(:,sys.jl1(i)));
    end
    
    % Specific treatment, where obj.uiq = unique(obj.iq) (see SystHBQ)
    Prods = (2*sys.vq)'.*( sys.Prod(Zr_compl(:,sys.jq),Zpmr_compl(:,sys.kq),'complex') );
    
    for i=1:numel(sys.uiq)
        Fpnl_fourier(:,sys.uiq(i))= Fpnl_fourier(:,sys.uiq(i)) + sum(Prods(:,(sys.iq==sys.uiq(i))),2);
    end
    
    %%% Does not seem to show better results because of the cost of unique function.
%     iq_update = obj.iq;
%     i_iq = obj.i_iq;
%     uiq = obj.uiq;
%     while ~isempty(uiq)
%         Fpnl_fourier(:,uiq)= Fpnl_fourier(:,uiq) + Prods(:,i_iq); % vectorial affectation.
%         Prods(:,uiq)=[];
%         iq_update(i_iq)=[];
%         [uiq,i_iq,~] = unique(iq_update);
%     end
    
    for i=1:numel(sys.id)
        Fpnl_fourier(:,sys.id(i))= Fpnl_fourier(:,sys.id(i)) + sys.vd(i)*(omegar*Zd_pmr(:,sys.jd(i))+omegapmr*Zd_r(:,sys.jd(i))) ;
    end
    
    for i=1:numel(sys.idd)
        Fpnl_fourier(:,sys.idd(i))= Fpnl_fourier(:,sys.idd(i)) + sys.vdd(i)*(omega2r*Zdd_pmr(:,sys.jdd(i))+omega2pmr*Zdd_r(:,sys.jdd(i))) ;
    end
    
    for i=1:numel(sys.id1)
        Fpnl_fourier(:,sys.id1(i))= Fpnl_fourier(:,sys.id1(i)) + sys.vd1(i)*(lambdaomegar*Zd_pmr(:,sys.jd1(i))+lambdaomegapmr*Zd_r(:,sys.jd1(i))) ;
    end
    
    % definition of omega2 and lambdaomega :
    aux_var = aux_var - [2*omegar*omegapmr; (lambdar*omegapmr + lambdapmr*omegar)];
    
    if ~isempty(sys.idl)
        for i=1:numel(sys.idq)
            Fpnl_fourier(:,sys.idq(i))= Fpnl_fourier(:,sys.idq(i)) + ...
                sys.vdq(i)*( sys.Prod(Zr_compl(:,sys.jdq(i)),sys.D(Zpmr_compl(:,sys.kdq(i)),'complex'),'complex') + ...
                sys.Prod(Zpmr_compl(:,sys.jdq(i)),sys.D(Zr_compl(:,sys.kdq(i)),'complex'),'complex'));
        end
        
        cond_ini(sys.idq)= cond_ini(sys.idq) + sys.vdq'.*((p-r)/p*Zr_t0(sys.jdq).*Zpmr_t0(sys.kdq) + r/p*Zpmr_t0(sys.jdq).*Zr_t0(sys.kdq));
        
    end
end


if mod(p,2) == 0
    Ur=get(U,'coefk',p/2);
    [Zr_tot,omegar,lambdar,omega2r,lambdaomegar] = sys.get_Ztot(Ur);
    
    % Vectorial computation of all the complex forms.
    Zr_compl = sys.real_to_compl(Zr_tot,'full');
    if ~isempty(sys.id) || ~isempty(sys.idd) || ~isempty(sys.id1) ; Zd_r = sys.D(Zr_tot); end
    if ~isempty(sys.idd) ; Zdd_r = sys.D(Zd_r); end
    
    Zr_t0 = sum(Zr_tot(1:H+1,:),1);
    
    % Loops on the lists.
    for i=1:numel(sys.iforce2)
        Fpnl_fourier(sys.hforce2(i),sys.iforce2(i))= Fpnl_fourier(sys.hforce2(i),sys.iforce2(i)) - sys.vforce2(i)*lambdar^2 ;
    end
    
    for i=1:numel(sys.ic2)
        Fpnl_fourier(1,sys.ic2(i))= Fpnl_fourier(1,sys.ic2(i)) + sys.vc2(i)*lambdar^2 ;
    end
    
    for i=1:numel(sys.il1)
        Fpnl_fourier(:,sys.il1(i))= Fpnl_fourier(:,sys.il1(i)) + (sys.vl1(i)*lambdar)*Zr_tot(:,sys.jl1(i));
    end
    
    % Specific treatment, where obj.uiq = unique(obj.iq) (see SystHBQ)
    Prods = sys.vq'.*( sys.Prod(Zr_compl(:,sys.jq),Zr_compl(:,sys.kq),'complex') );
    
    for i=1:numel(sys.uiq)
        Fpnl_fourier(:,sys.uiq(i))= Fpnl_fourier(:,sys.uiq(i)) + sum(Prods(:,(sys.iq==sys.uiq(i))),2);
    end
    
    %%% Does not seem to show better results because of the cost of unique function.
%     iq_update = obj.iq;
%     i_iq = obj.i_iq;
%     uiq = obj.uiq;
%     while ~isempty(uiq)
%         Fpnl_fourier(:,uiq)= Fpnl_fourier(:,uiq) + Prods(:,i_iq); % vectorial affectation.
%         Prods(:,uiq)=[];
%         iq_update(i_iq)=[];
%         [uiq,i_iq,~] = unique(iq_update);
%     end
    
    for i=1:numel(sys.id)
        Fpnl_fourier(:,sys.id(i))= Fpnl_fourier(:,sys.id(i)) + (sys.vd(i)*omegar)*Zd_r(:,sys.jd(i)) ;
    end
    
    for i=1:numel(sys.idd)
        Fpnl_fourier(:,sys.idd(i))= Fpnl_fourier(:,sys.idd(i)) + (sys.vdd(i)*omega2r)*Zdd_r(:,sys.jdd(i)) ;
    end
    
    for i=1:numel(sys.id1)
        Fpnl_fourier(:,sys.id1(i))= Fpnl_fourier(:,sys.id1(i)) + (sys.vd1(i)*lambdaomegar)*Zd_r(:,sys.jd1(i)) ;
    end
    
    % definition of omega2 and lambdaomega :
    aux_var = aux_var - [omegar^2; lambdar*omegar];
    
    if ~isempty(sys.idl)
        for i=1:numel(sys.idq)
            Fpnl_fourier(:,sys.idq(i))= Fpnl_fourier(:,sys.idq(i)) + ...
                sys.vdq(i)*( sys.Prod(Zr_compl(:,sys.jdq(i)),sys.D(Zr_compl(:,sys.kdq(i)),'complex'),'complex'));
        end
        
        cond_ini(sys.idq)= cond_ini(sys.idq) + sys.vdq'.*(1/2*Zr_t0(sys.jdq).*Zr_t0(sys.kdq));
    end
    
end
    
% transform the matrix Fpnl to a column vector and add 0 for the phase
% equation, omega2,.. eqautions
Fpnl_fourier(1,sys.idq) = cond_ini(sys.idq);
Fpnl_fourier = reshape(Fpnl_fourier,nz_tot*DHp1,1);

Fpnl = [Fpnl_fourier(1:sys.neq-1);Fpnl_phase];
Fpnl_aux = [Fpnl_fourier(sys.neq:end);aux_var];
Fpnl_tot = [Fpnl;Fpnl_aux];
