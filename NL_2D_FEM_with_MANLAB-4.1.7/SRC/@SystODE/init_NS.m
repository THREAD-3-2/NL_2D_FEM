function U0 = init_NS(sys,Section)

%Utemp = evalseries(Section.U0,Apoint,sys.order);
Utemp = Section.Ustab;
Hold = ((numel(Utemp)-4)/sys.nz_tot-1)/2;
H = sys.H;

omega1 = Utemp(sys.nz*(2*Hold+1)+1);
lambda = Utemp(sys.nz*(2*Hold+1)+2);
Ztemp = reshape(Utemp([1:sys.nz*(2*Hold+1) sys.nz*(2*Hold+1)+3:sys.nz_tot*(2*Hold+1)+2]),2*Hold+1,sys.nz_tot);

% truncation of Ztemp to H<=Hold harmonics
ind_kept = [1:H+1 Hold+2:Hold+H+1];
Zper = Ztemp(ind_kept,:);

coef_per_var = (2*(2*H+2)*H+1);
Eigen = Section.Eigen;

[~,ind] = sort(abs(real(Eigen.values)));
Eigen.values = Eigen.values(ind);

Eigen.vectors = Eigen.vectors(:,ind);
omega2 = imag(Eigen.values(1));
if omega2<0
    Eigen.vectors(:,[1 2]) = Eigen.vectors(:,[2 1]);
    Eigen.values([1 2]) = Eigen.values([2 1]);
    omega2 = abs(omega2);
end

Pert_w2 = reshape(Eigen.vectors(:,1),2*Hold+1,sys.nz_tot);
%%% Second phase equation :
% Attention ce n'est pas tout à fait ça.. !
%Pert_w2 = conj(Pert_w2(sys.zi_phase2,1))/abs(Pert_w2(sys.zi_phase2,1)) * Pert_w2(ind_kept,:);
%%% Normalization
valno = 5e-2*norm(Zper);
val_normalize = valno / max(abs(reshape(Pert_w2,numel(Pert_w2),1)));
Pert_w2 = Pert_w2*val_normalize;

DHp1 = 2*H+1;
ind_change_cos = reshape([DHp1 DHp1+2] + (0:1:(H-1))'*DHp1,DHp1-1,1);
ind_change_sin = ind_change_cos+(coef_per_var-1)/2;
ind_cos_per = DHp1+1 + (0:1:(H-1))*DHp1;
ind_sin_per = ind_cos_per+(coef_per_var-1)/2;

Pck = Pert_w2(2:H+1,:);
Psk = Pert_w2(Hold+2:Hold+H+1,:);

ZNS = randn(coef_per_var,sys.nz_tot)*(norm(Zper)*1e-7);
ZNS([1 ind_cos_per ind_sin_per],:) = Zper;
ZNS(2,:) = real(Pert_w2(1,:));
ZNS((coef_per_var+1)/2+1,:) = -imag(Pert_w2(1,:));
ZNS(ind_change_cos,:) = 0.5*[real(Pck)-imag(Psk) ; real(Pck)+imag(Psk)];
ZNS(ind_change_sin,:) = 0.5*[imag(Pck)+real(Psk) ; imag(Pck)-real(Psk)];

lambda = (1+1e-2)*lambda;
if strcmp(sys.subtype,'forced')
    omega1 = lambda;
end
    
U0 = sys.init_U0(ZNS,omega1,omega2,lambda);

end