function [Utf,Vtf,Jacobian] = tangentvector(sys,Jacobian)
% Computation of the tangent vector at U0 using the Jacobian matrix dRdU at U0
% dRdU is augmented with a random vector

% We assume that lambda is not a limit-point. In theory there could be
% issues but the ending point of a continuation branch is never close
% enough to a limit-point

% It returns the tangent vector before normalization : to
% improve efficiency of obj.ANMSeries.

%% Test if there are auxiliary variables.
if ~isempty(Jacobian.dRauxdUaux)
    [Jacobian.LK_aux,Jacobian.UK_aux,Jacobian.pK_aux,Jacobian.qK_aux] = lu(Jacobian.dRauxdUaux,'vector');
    
%    Juka = Jacobian.UK_aux; Jlka=Jacobian.LK_aux; Jdradu = Jacobian.dRauxdU(Jacobian.pK_aux,:);
    
   Jacobian.Kfull = Jacobian.dRdU - ...
       Jacobian.dRdUaux(:,Jacobian.qK_aux)*(Jacobian.UK_aux\(Jacobian.LK_aux\Jacobian.dRauxdU(Jacobian.pK_aux,:)));

%     Jacobian.Kfull = sparse(Jacobian.dRdU - ...
%         Jacobian.dRdUaux(:,Jacobian.qK_aux)*(pinv(full(Jacobian.UK_aux))*(Jacobian.LK_aux\Jacobian.dRauxdU(Jacobian.pK_aux,:))));

else
    Jacobian.Kfull = Jacobian.dRdU;
end
%% Default column we put in the right-hand-side
Jacobian.missingcolumn = Jacobian.Kfull(:,end);
Jacobian.K = Jacobian.Kfull(:,1:end-1);
% Fast lu factorization
[Jacobian.LK,Jacobian.UK,Jacobian.pK,Jacobian.qK] = lu(Jacobian.K,'vector');
%% To avoid issues with vertical branches
% The column sent to the right-hand-side changes if
% an eigenvalue of K is smaller than tol.
[valmin,imin] = min(abs(diag(Jacobian.UK)));
[valmax,imax] = max(abs(diag(Jacobian.UK)));
val = valmin/valmax;
tol = 1e-12;
if val < tol   
    if valmax > 1/valmin
        i = imax;
    else
        i = imin;
    end
    % The column corresponding to the smallest eigenvalue is sent to the
    % right-hand side.
    Jacobian.index_exchange = Jacobian.qK(i);
    Jacobian.list_val = [1:Jacobian.index_exchange-1 Jacobian.index_exchange+1:sys.neq+1];
    Jacobian.missingcolumn = Jacobian.Kfull(:,Jacobian.index_exchange);
    Jacobian.K = Jacobian.Kfull(:,Jacobian.list_val);
    % Fast lu factorization
    [Jacobian.LK,Jacobian.UK,Jacobian.pK,Jacobian.qK] = lu(Jacobian.K,'vector');
    %disp('There is a vertical tangent. Computation may be only approximations.');
    Jacobian.vertical_tangent = 1;
else
    Jacobian.index_exchange = sys.neq+1;
    Jacobian.list_val = 1:sys.neq;
    
    Jacobian.vertical_tangent = 0;
end
[Vt,Vt_aux] = sys.Condensation(Jacobian,sparse(sys.neq,1),sparse(sys.neq_aux,1),1);
Vtf = [Vt;Vt_aux];

Utf = Vtf *  1/sqrt((sys.arclengthdef.*Vtf)'*Vtf);


end