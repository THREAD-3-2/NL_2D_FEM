function [eig_values,eig_vectors] = Stability(sys,Utot,Jacobian)
% function [eig_values,eig_vectors] = Stability(obj,Utot,Jacobian)
% Computes the stability of a solution Utot with the Jacobian
% matrix K using the appropriate method.

switch sys.type
    case 'Equ'
        [eig_values,eig_vectors]=sys.Stability_Equ(Utot,Jacobian);
    case 'HBM'
        [eig_values,eig_vectors]=sys.Stability_HBM(Utot,Jacobian);
    case 'QPHBM'
        [eig_values,eig_vectors]=sys.Stability_QPHBM(Utot,Jacobian);
%         warndlg('Stability is not available for Quasi-periodic solutions.');
         eig_values = NaN;
         eig_vectors = NaN;
end
