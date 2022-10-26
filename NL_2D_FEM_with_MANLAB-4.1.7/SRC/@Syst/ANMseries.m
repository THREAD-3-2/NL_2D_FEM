function [U0_tot,Amax, BifData, StabData] = ANMseries(sys,U0tot,Ut0)
% Compute ANM series for sys(U0val)
% Outputs: U0_tot (taylor type), Amax (domaine of utility),
% BifData is a structure variable used to store bifurcation point data
% when one is found  properties : Ubif,Utf, Utb,alpha,status
% StabData is a structure variable used to store stability data
% when a change of stability is found  properties : Ustab, Astab, status

global Ck

%% Tangent vector - branch orientation

Jacobian = sys.Jacobian(U0tot);% Block decomposition of the Jacobian at point U0tot

[U1,Vt_tot,Jacobian] = sys.tangentvector(Jacobian);     % tangent vector at U0tot

if (U1'*Ut0) < 0, U1 = - U1; end  % set U1 in the direction of Ut0

%% series computation : Zero and First order
U0_tot=Taylor(get(sys,'order'),U0tot);
U0_tot=set(U0_tot,'coef1',U1);

%% series computation : Orders p>=2
Rhs0_tot = U1/(U1'*Vt_tot);

warning('off','MATLAB:nearlySingularMatrix');
warning('off','MATLAB:singularMatrix');
for p=2:sys.order
    Ck=p;
    
    % Fast computation using the condensation function :
    Fpnl_tot = sys.Fpnl(U0_tot,p);
    
    [Rhs1,Rhs1_aux] = Condensation(sys,Jacobian,-Fpnl_tot(1:sys.neq),-Fpnl_tot(sys.neq+1:end),0);
    Rhs1_tot = [Rhs1;Rhs1_aux];
    
    Usol_missing = - Rhs0_tot'*Rhs1_tot;
    Usol_tot = Rhs1_tot + Usol_missing*Vt_tot;
    
    U0_tot=set(U0_tot,'coefk',Usol_tot,p);
end

%% Detection and Extraction of emerging geometric series

if sys.BifDetection == 1
    [Uclean, BifData.Uscale, BifData.alpha, BifData.status] = sys.GeomSerie(U0_tot);
    
    switch BifData.status
        case 'nothing'
            
        case 'simplebif'
            U0_tot=Uclean;
            Fpnl_tot=sys.Fpnl(U0_tot,p);
            [BifData.Ubif, BifData.Utf, BifData.Utb]=sys.locatebif(Uclean, BifData.Uscale, BifData.alpha);
    end
else
    BifData.Uscale=[];
    BifData.alpha=0;
    BifData.status = 'nothing';
end

%% domaine of utility Amax
if (norm(Fpnl_tot)==0)
    Amax = sys.Amax_max;
    disp(['ANMSeries:Infinite radius, Amax set to ' num2str(Amax)]);
  else
    Amax = (sys.ANMthreshold/norm(Fpnl_tot))^(1/(sys.order));
    Amax = min(Amax,sys.Amax_max);
end

%% Stability of the solution
if sys.StabilityCheck == 1
    StabData = sys.StabilityComputation(U0_tot,Amax);
else
    StabData.Eigen_init = 0;
    StabData.Eigen_end = 0;
    StabData.Uchange = 0;
    StabData.Achange = 0;
    StabData.Eigen.type = 'nothing';
    StabData.status = {'stable','stable'};
end

end

