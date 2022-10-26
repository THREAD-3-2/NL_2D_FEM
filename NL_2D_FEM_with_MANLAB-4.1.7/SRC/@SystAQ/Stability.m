function [eig_val,eig_vec] = Stability(sys,Utot,Jacobian)
% function [status,realpos,eig_val,eig_vec] = Stability(obj,Utot,Jacobian)
% An additional argument can be specify as an initial guess of the
% eigenvectors and eigenvalues.

%% Linear stability
K = Jacobian.K;

[eig_vec,eig_val] = eig(full(K));
eig_val = diag(eig_val);

[~,ind] = sort(abs(real(eig_val)));
eig_val = eig_val(ind);
eig_vec = eig_vec(:,ind);


if ~isempty(Jacobian.dRauxdUaux)
    
    Kaux = Jacobian.dRauxdUaux;
    
    %%% Without lambda
    eig_vec_aux = - Kaux\(Jacobian.dRauxdU(:,1:end-1)*eig_vec);
    eig_vec = [eig_vec ; eig_vec_aux];
end



end


