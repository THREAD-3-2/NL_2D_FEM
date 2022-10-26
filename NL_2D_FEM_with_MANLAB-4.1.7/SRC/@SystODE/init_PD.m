function [Unew] = init_PD(sys,Section)
% function [Unew] = init_PD(sys,Section)
% Initialize at a period doubling point by adding zero in the 
% (2k+1) * (omega/2) terms of the Fourier series.

Uold = Section.Ustab;
Hold = ((size(Uold,1)-4)/sys.nz_tot-1)/2;

neq_old = sys.nz*(2*Hold+1);
neq_aux_old = sys.nz_aux*(2*Hold+1);
Zold = reshape(Uold([1:neq_old neq_old+3:neq_old+3+neq_aux_old-1]) , 2*Hold+1 , sys.nz_tot); 
omega = Uold(neq_old+1)/2;
lambda = Uold(neq_old+2);

H = sys.H;
DHp1 = 2*H+1;

Znew = zeros(DHp1,sys.nz_tot);
Znew([1 (3:2:min(2*Hold+2,H+1)) (H+3:2:min(2*H+1,H+2+2*Hold))],:) = Zold([(1:min(H/2+1,Hold+1)) (Hold+2:Hold+1+min(H/2,Hold))],:);

Unew = sys.init_U0(Znew,omega,lambda);

end


% % % Try with Hill eigenvectors. Does not work.
% 
% Utemp = Section.Ustab;
% Hold = ((numel(Utemp)-4)/sys.nz_tot-1)/2;
% H = sys.H;
% 
% omega = Utemp(sys.nz*(2*Hold+1)+1);
% lambda = Utemp(sys.nz*(2*Hold+1)+2);
% Ztemp = reshape(Utemp([1:sys.nz*(2*Hold+1) sys.nz*(2*Hold+1)+3:sys.nz_tot*(2*Hold+1)+2]),2*Hold+1,sys.nz_tot);
% 
% Ztemp_compl = sys.real_to_compl(Ztemp);
% 
% % truncation of Ztemp to H<=Hold harmonics
% Ztemp_compl = Ztemp_compl(1:floor(H/2)+1,:);
% 
% Eigen = Section.Eigen;
% [~,ind] = sort(abs(real(Eigen.values)));
% Eigen.vectors = Eigen.vectors(:,ind);
% 
% Zpert = reshape(Eigen.vectors(:,1),2*Hold+1,sys.nz_tot);
% valnormalization = 1e-2*norm(Ztemp_compl);
% 
% Zpert = valnormalization*Zpert;
% Zpert_compl = sys.real_to_compl(Zpert);
% Zpert_compl = Zpert_compl(1:floor(H/2),:);
% 
% Zpd_compl = zeros(H+1,sys.nz_tot);
% Zpd_compl(1:2:H+1,:) = Ztemp_compl;
% Zpd_compl(2:2:H+1,:) = Zpert_compl;
% 
% phase = conj(Zpd_compl(2,sys.zi_phase1))/abs(Zpd_compl(2,sys.zi_phase1));
% for ii=1:H
%     Zpd_compl(1+ii,:) = Zpd_compl(1+ii,:)*phase^ii;
% end
% 
% Zpd = sys.compl_to_real(Zpd_compl);
% 
% if strcmp(sys.subtype,'forced'); lambda = omega/2; disp('Be careful with the writing of the forcing'); end
% 
% U0 = sys.init_U0(Zpd,omega/2,lambda); 
% 
% end
% 
