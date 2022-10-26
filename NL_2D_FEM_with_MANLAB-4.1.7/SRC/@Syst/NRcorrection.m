function [Uf] = NRcorrection(sys,Uf)
% performs Newton-Raphson correction

iter = 0;

NRthreshold=get(sys,'NRthreshold'); NRitemax=get(sys,'NRitemax');
NRmethod = get(sys,'NRmethod');
Usave = Uf;
Rnew_f=sys.R(sys,Uf) + sys.pertvect*sys.PerturbationSize;
normR =norm(Rnew_f);
best_normR = normR;
if (normR > NRthreshold)
    disp(['N-R: Before correction ||R(Uj)||=' num2str(normR); ] );
end

while (normR > NRthreshold) && (iter < NRitemax)
    iter = iter + 1;
    Jacobian = sys.Jacobian(Uf);
    
    % Check if Utot is a limit-point and update the structure Jacobian.
    [Utf,~,Jacobian] = sys.tangentvector(Jacobian);
    
    if NRmethod == 1 % The parameter is fixed, there are no additional inversion of the Jacobian matrix needed.
        % Corrections thanks to the condensation function
        [Ucor,Ucor_aux] = sys.Condensation(Jacobian,Rnew_f(1:sys.neq),Rnew_f(sys.neq+1:end),0);
        Ucor_f = [Ucor;Ucor_aux];
        
        Uf=Uf-Ucor_f;
                
    elseif (NRmethod == 2 || NRmethod == 0)% The parameter is variable, an additional inversion of the Jacobian matrix is needed.
        
        Jacobian.Kfull = [Jacobian.Kfull ; Utf(1:sys.neq+1)'];
        if sys.neq_aux > 0
            Jacobian.dRdUaux = [Jacobian.dRdUaux ; Utf(sys.neq+2:end)'];
            
            B = [Rnew_f(1:sys.neq);0];
            Baux = Rnew_f(sys.neq+1:end);
            % Right-hand-side / Condensation formulaes
            Rhs = B - Jacobian.dRdUaux(:,Jacobian.qK_aux)*(Jacobian.UK_aux\(Jacobian.LK_aux\Baux(Jacobian.pK_aux)));
            
            % Computation of the main variables, then the auxiliary variables :
            Um = Jacobian.Kfull\Rhs;
            Ua = zeros(sys.neq_aux,1);
            Ua(Jacobian.qK_aux) = Jacobian.UK_aux\(Jacobian.LK_aux\(Baux(Jacobian.pK_aux) - Jacobian.dRauxdU(Jacobian.pK_aux,:)*Um));
            
            Ucor_f = [Um;Ua];
        else
            Ucor_f = Jacobian.Kfull\[Rnew_f;0];
        end
        
        Uf=Uf-Ucor_f;
    end
    
    Rnew_f=sys.R(sys,Uf) + sys.pertvect*sys.PerturbationSize;
    normR = norm(Rnew_f);
    if normR < best_normR
        best_normR = normR;
        Usave = Uf;
    end
    disp(['N-R: After correction ' num2str(iter) ' : ||R(Uj)||=' num2str(normR); ] );
end

if isnan(normR)
    warndlg('Correction failed: NaN values. Back to the best point.','Correction','on');
    Uf = Usave;
end

if iter == NRitemax
    disp('Correction failed: Max number of iterations is reached. Back to the best point.');
    Uf = Usave;
end
