function [U, Amax,BifData, StabData]=ANMserieBif(sys,Ubif,Ut1,Ut2)
%Compute taylor series at Ubif for the branch along Ut1 tangent

global Ck
Ck=sys.order;

Norder = get(sys,'order');

Jacobian=sys.Jacobian(Ubif);

% Transpose of the Jacobian matrix :
TJacobian.dRdU = Jacobian.dRdU';
TJacobian.dRauxdUaux = Jacobian.dRauxdUaux';
if any(Jacobian.dRauxdUaux)
    TJacobian.dRdUaux = Jacobian.dRauxdU';
    TJacobian.dRauxdU = Jacobian.dRdUaux';
    
    [TJacobian.LK_aux,TJacobian.UK_aux,TJacobian.pK_aux,TJacobian.qK_aux] = lu(TJacobian.dRauxdUaux,'vector');
    TJacobian.Kfull = TJacobian.dRdU - ... 
        TJacobian.dRdUaux(:,TJacobian.qK_aux)*(TJacobian.UK_aux\(TJacobian.LK_aux\TJacobian.dRauxdU(TJacobian.pK_aux,:)));

    TJacobian.dRdUaux = TJacobian.dRdUaux(1:end-1,:);
else
    TJacobian.Kfull = TJacobian.dRdU;
end

TJacobian.K = TJacobian.Kfull(1:end-1,:);
[TJacobian.LK,TJacobian.UK,TJacobian.pK,TJacobian.qK] = lu(TJacobian.K,'vector');

% left nul vector phi computed with two inverse iterations using condensation
[phi,phi_aux] = sys.Condensation(TJacobian,sparse(ones(sys.neq,1)),sparse(ones(sys.neq_aux,1)));
phi_tot = [phi;phi_aux];
phi_tot = phi_tot/norm(phi_tot);
[phi,phi_aux] = sys.Condensation(TJacobian,phi_tot(1:sys.neq),phi_tot(sys.neq+1:end));
phi_tot = [phi;phi_aux];
phi_tot = phi_tot'/norm(phi_tot);

phi = phi_tot(1:sys.neq);
phi_aux = phi_tot(sys.neq+1:end);

%Nul right vector of the matrix [dRdU; Ut1']
Umode = Ut2 - ( (Ut1'*Ut2)/(Ut1'*Ut1)  ) *Ut1;

% Augmented jacobian matrix :
Jacobian.dRdU = [ Jacobian.dRdU phi'; [Ut1(1:sys.neq+1) ; 0]' ;[Umode(1:sys.neq+1); 0]'];
if ~isempty(phi_aux)
    Jacobian.dRauxdU = [Jacobian.dRauxdU phi_aux'];
    Jacobian.dRdUaux = [Jacobian.dRdUaux ; Ut1(sys.neq+2:end)' ; Umode(sys.neq+2:end)' ];

    [Jacobian.LK_aux,Jacobian.UK_aux,Jacobian.pK_aux,Jacobian.qK_aux] = lu(Jacobian.dRauxdUaux,'vector');
    Jacobian.K = Jacobian.dRdU - ... 
        Jacobian.dRdUaux(:,Jacobian.qK_aux)*(Jacobian.UK_aux\(Jacobian.LK_aux\Jacobian.dRauxdU(Jacobian.pK_aux,:)));
else
    Jacobian.K = Jacobian.dRdU;
end
[Jacobian.LK,Jacobian.UK,Jacobian.pK,Jacobian.qK] = lu(Jacobian.K,'vector');

% Initialization of Taylor object
U = Taylor(Norder+1,Ubif);

% Order 1
U = set(U, 'coef1', Ut1);

% prepare some constant
U=set(U,'coefk',Umode, 2);
U=set(U,'coefk',zeros(size(Ubif)), 3);
Ck=3;
F12_tot = sys.Fpnl(U,3);
ProjF12= phi_tot*F12_tot;


% Second Order
% ------------
%rhs for order 2

U=set(U,'coefk',zeros(size(Ubif)), 2);
Ck=2;
Fpnl_tot    =  sys.Fpnl(U,2);
[Usol,Usol_aux] = sys.Condensation(Jacobian,-[Fpnl_tot(1:sys.neq);0;0],-Fpnl_tot(sys.neq+1:end));

U2star= [Usol(1:sys.neq+1);Usol_aux];
U=set(U,'coefk',U2star, 2);

% Compute alpha, finish U2 and prepare rhs for order 3

U=set(U,'coefk',zeros(size(Ubif)), 3);
Ck=3;
F3nlE_tot = sys.Fpnl(U,3);
alpha =  - (phi_tot*F3nlE_tot)/ProjF12;
U2    =  U2star + alpha * Umode;
Fknl_tot  = F3nlE_tot + alpha * F12_tot;

U =set(U,'coefk',U2,2); %

%%
% Next Orders
%------------

for p=3:Norder
    
    % Use the condensation formulaes :
    [Usol,Usol_aux] = sys.Condensation(Jacobian,-[Fknl_tot(1:sys.neq);0;0],-Fknl_tot(sys.neq+1:end));
    Upstar= [Usol(1:sys.neq+1);Usol_aux];
    U =set(U,'coefk',Upstar,p);
    
    % Compute alpha, finish Up and prepare rhs for order p+1
    U =set(U,'coefk',zeros(size(Ubif)),p+1);
    Ck=p+1;
    FknlE_tot = sys.Fpnl(U,p+1);
    alpha =  - (phi_tot*FknlE_tot)/ProjF12;
    Up    =  Upstar + alpha * Umode;
    Fknl_tot  = FknlE_tot + alpha * F12_tot;
    U =set(U,'coefk',Up,p); %
end

Ck = Norder; % Go back to the order of the series.

BifData.status='nothing';

%% domaine of utility Amax
if (norm(Fknl_tot)==0)
    Amax = sys.Amax_max;
    disp(['ANMSeries:Infinite radius, Amax set to ' num2str(Amax)]);
  else
    Amax = (sys.ANMthreshold/norm(Fknl_tot))^(1/(sys.order+1));
    Amax = min(Amax,sys.Amax_max);
end

%% Stability of the solution
if sys.StabilityCheck == 1
    StabData = sys.StabilityComputation(U,Amax);
else
    StabData.Eigen_init = 0;
    StabData.Eigen_end = 0;
    StabData.Uchange = 0;
    StabData.Achange = 0;
    StabData.Eigen.type = 'nothing';
    StabData.status = {'stable','stable'};
end

