function [eig_val,eig_vec] = Stability_Equ(sys,U,Jacobian)
% function [eig_val,eig_vec] = Stability_Equ(sys,Ut,Jacobian)
% Computes the eigenvalue and the eigenvectors of 

%%% Generalization with mass matrix

Massfull = sparse(sys.id,sys.jd,-sys.vd,sys.nz_tot,sys.nz_tot);
Mass = Massfull(1:sys.nz,1:sys.nz);
Massaux = Massfull(sys.nz+1:end,1:sys.nz);

Kaux = Jacobian.dRauxdUaux;
B = Jacobian.dRdUaux;

MJ = Mass - B*(Kaux\Massaux);

%%% Computation of eigenvalues and eigenvectors
% [eig_vec,eig_val] = eig(full(Jacobian.K)); % Version without mass matrix
[eig_vec,eig_val] = eig(full(Jacobian.K),full(MJ));
[eig_val] = diag(eig_val);

%%% Detection of infinite eigenvalues
if sum(isinf(eig_val))>0
    disp('SystODE/Stability_Equ.m : Some eigenvalues seems infinite. Stability results might be wrong.');
end

%%% Sorting
[~,ind] = sort(abs(real(eig_val)));
eig_val = eig_val(ind);
eig_vec = eig_vec(:,ind);

%%% Auxiliary eigenvectors
if ~isempty(Jacobian.dRauxdUaux)
    
    Kaux = Jacobian.dRauxdUaux;
    C = Jacobian.dRauxdU;
    
    %%% Without lambda
    eig_vec_aux = - Kaux\(C(:,1:end-1)*eig_vec);
    eig_vec = [eig_vec ; eig_vec_aux];
end