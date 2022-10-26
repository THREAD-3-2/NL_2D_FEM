function [floquet_exponent,hill_eigenvectors] = Stability_HBM(sys,Utot,Jacobian)
%function [floquet_exponent,hill_eigenvectors] = Stability(sys,Utot,Jacobian)
% Computes the stability of a periodic solution Utot with the Jacobian
% matrix K using Hill's purely frequency based method

K = Jacobian.K;

% The Hill matrix does not include the angfreq contribution
Hill_mat = K(1:end-1,1:end-1);

% Generalization with mass matrix
M = Jacobian.M;
Maux = Jacobian.Maux;

MassHill = M - Jacobian.dRdUaux(:,Jacobian.qK_aux)*(Jacobian.UK_aux\(Jacobian.LK_aux\Maux(Jacobian.pK_aux,:)));
MassHill = MassHill(1:end-1,1:end-2);

% Diagonalization of Hill's matrix
% [V,D] = eig(full(Hill_mat)); % Version without mass matrix
[V,D] = eig(full(Hill_mat),full(MassHill));
eig_values = diag(D);

%%% Detection of infinite eigenvalues
if sum(isinf(eig_values))>0
    disp('SystODE/Stability_HBM.m : Some eigenvalues are infinite. Stability results might be wrong.');
end

% Computation of the complex eigenvectors
H = sys.H;
DHp1 = 2*H+1;
nz = sys.nz;

ind_cos = 2:H+1;
ind_costot = repmat(ind_cos,nz,1)+repmat((0:DHp1:(nz-1)*DHp1)',1,H);
ind_sin = H+2:DHp1;
ind_sintot = repmat(ind_sin,nz,1)+repmat((0:DHp1:(nz-1)*DHp1)',1,H);
Vcst = V(1:DHp1:end,:);
Vcos = V(ind_costot,:);
Vsin = V(ind_sintot,:);
Vexp_pos = (Vcos - 1i*Vsin)/2;
Vexp_neg = (Vcos + 1i*Vsin)/2;
Vcomp = [Vexp_neg(end:-1:1,:);Vcst;Vexp_pos];

% Sorting of the eigenvector to keep the most converged
XX = reshape(ones(nz,1)*(-H:H),nz*(2*H+1),1); % used for sorting
Wmean = (sum(repmat(XX,1,nz*(2*H+1)).*abs(Vcomp),1)./sum(abs(Vcomp),1));
[val,ind1] = sort(Wmean);
ind2 = find((val < 0.5+sys.StabTol) & (val > -0.5+sys.StabTol));            %%% BUG


[~,ind3] = sort(abs(real(eig_values(ind1(ind2)))));
floquet_exponent = eig_values(ind1(ind2(ind3)));
hill_eigenvectors = V(:,ind1(ind2(ind3))); % For later use

% Test of convergence, nz values should be between -0.5 and 0.5 
if (length(ind2) < nz/10 || length(ind2) > nz)
    disp('SystODE/Stability_HBM.m : Floquet Exponents not converged. Stability results may be wrong.');
end

% If the system is autonomous, the null floquet exponent is removed. 
if strcmp(sys.subtype,'autonomous') %% There should be one floquet exponent which is null.
    [~,ind_suppr] = min(abs(floquet_exponent));
    hill_eigenvectors(:,ind_suppr) = [];
    floquet_exponent(ind_suppr) = [];
end

%% Auxiliary eigenvectors ; to compute emerging quasi-periodic solution at NS bifurcation.
if ~isempty(Jacobian.dRauxdUaux) 
    Kaux = Jacobian.dRauxdUaux;
    C = Jacobian.dRauxdU;
    
    %%% Without omega,lambda | omega^2, lambda*omega
    hill_eigvec_aux = - Kaux(1:end-2,1:end-2)\(C(1:end-2,1:end-2)*hill_eigenvectors);
    hill_eigenvectors = [hill_eigenvectors ; hill_eigvec_aux];
    
    %%% The auxiliary eigenvectors may need something else if time
    %%% derivative of the main variables appear in auxiliary equations
end

% figure(7)
% plot(floquet_exponent,'o');hold on;
% line([0 0],[min(imag(floquet_exponent)) max(imag(floquet_exponent))],'color','k');


end


