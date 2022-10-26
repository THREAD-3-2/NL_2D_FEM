function [Ubif, Utf, Utb]=locatebif(sys, Uclean,Uscale,Alpha)
% locate a simple bifurcation point and its two tangents
% outputs : Ubif (bif point), Utf and Utb (the two tangents
% at Ubif along fondamental path and bifurcated path)
% inputs : Alpha (common ratio), Uscale (scale factor)
global Ck
Ck=sys.order;

Norder=get(sys, 'order');

Ubif=evalseries(Uclean,Alpha,Norder);

Jacobian=sys.Jacobian(Ubif);

% Transpose of the Jacobian matrix :
TJacobian.dRtotdUtot = Jacobian.dRtotdUtot';
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

%Utf: tangent of the fondamental path.
Utf=evalderiv(Uclean,Alpha,Norder);

%normalisation of Utf and Uscale
Utf=Utf/norm(Utf);
Uscale=Uscale/norm(Uscale);

%Computation of the second tangent at the biffurcation point
% a) Computation of ddRddU_Utb(Uscale,Uscale)
Utemp=Taylor(Norder,Ubif);
Utemp1=set(Utemp,'coef1',Uscale); % Initialization of first derivative U1 to Utf
Utemp1=set(Utemp1,'coefk',zeros(size(Utf)),2);
ddRUdd_UU=sys.Fpnl(Utemp1,2);

% b) computation of ddRddU(Utf-Uscale) and ddRddU(Utf+Uscale)
Vtemp=set(Utemp1, 'coef1', Utf + Uscale);
Wtemp=set(Utemp1, 'coef1', Utf - Uscale);

ddRUdd_UV=(sys.Fpnl(Vtemp,2)-sys.Fpnl(Wtemp,2))./2;

Utb=Utf-(phi_tot*ddRUdd_UV)/(phi_tot*ddRUdd_UU)*Uscale;
Utb=Utb/norm(Utb);

% check accuracy of the results
Rbif_tot=sys.R(sys,Ubif);
Rtf=Jacobian.dRtotdUtot*Utf;
Rtb=Jacobian.dRtotdUtot*Utb;

%phi'*Q(obj,Utf,Utf)
%Utemp1=set(Utemp1,'coef1',Utf);
%eqbiff=phi_tot*obj.Fpnl(Utemp1,2);

%phi'*Q(obj,Utb,Utb)
%Utemp1=set(Utemp1,'coef1',Utb);
%eqbifb=phi_tot*obj.Fpnl(Utemp1,2);

disp(['Bif info: ||R(Ubif)||=' num2str(norm(Rbif_tot))]);
%disp(['Left-nul vector: ||TdRdU(TPhi)||=' num2str(norm(TJacobian.dRtotdUtot*phi_tot')) ' ||phi||=' num2str(norm(phi_tot))]);
disp(['Tangent fond: ||dRdU(Utf)||=' num2str(norm(Rtf)) ' ||Utf||=' num2str(norm(Utf))]);
disp(['Tangent bif : ||dRdU(Utb)||=' num2str(norm(Rtb)) ' ||Utb||=' num2str(norm(Utb))]);
%disp(['Bif info: TPhi.Q(Utf,Utf)=' num2str(eqbiff) ' TPhi.Q(Utb, Utb)=' num2str(eqbifb)]);
disp(' ');

