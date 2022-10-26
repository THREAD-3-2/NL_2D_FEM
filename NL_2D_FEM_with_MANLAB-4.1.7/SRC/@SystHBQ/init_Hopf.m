function U0 = init_Hopf(sys,Section)

Eigen = Section.Eigen;
[~,ind] = sort(abs(real(Eigen.values)));
Eigen.values = Eigen.values(ind);
Eigen.vectors = Eigen.vectors(:,ind);

Ucst = Section.Ustab;

omega = imag(Eigen.values(1));
if omega<0
    Eigen.vectors(:,[1 2]) = Eigen.vectors(:,[2 1]);
    Eigen.values([1 2]) = Eigen.values([2 1]);
    omega = abs(omega);
end

% lambda was in neq+1 = nz+1 position
lambda = Ucst(sys.nz+1);
% lambda is removed from the stationary part.
Ucst = Ucst([1:sys.nz sys.nz+2:sys.nz_tot+1]);

% Constant values initialized to the stationary solution
Ucompl_H1 = [Ucst ; conj(Eigen.vectors(sys.zi_phase1,1))/abs(Eigen.vectors(sys.zi_phase1,1))*Eigen.vectors(:,1)];
Zreal_H1 = sys.compl_to_real(transpose(reshape(Ucompl_H1,sys.nz_tot,2)));

% Normalization of the arising periodic orbit
valno = 5e-2;
val_normalize = valno/ max(abs([Zreal_H1(2,:) Zreal_H1(3,:)]));

% Arising periodic orbit approximation
ZHopf = randn(2*sys.H+1,sys.nz_tot)*1e-5;
ZHopf([1 2 2+sys.H],:) = Zreal_H1 * val_normalize;

U0 = sys.init_U0(ZHopf,omega,lambda);

end