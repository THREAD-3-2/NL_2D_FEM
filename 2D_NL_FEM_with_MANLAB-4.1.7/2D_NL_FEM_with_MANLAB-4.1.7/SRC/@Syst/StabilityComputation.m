function StabData = StabilityComputation(sys,U0,Amax)
% function StabData = StabilityComputation(obj,U0,Amax)
% This function computes the point of changement of stability if there is a
% change of stability. This function assumes that there is only one change
% of stability of the solution, otherwise, it can give erroneous results.
%
% The exact position is computed using polynomial fits on the real parts of
% the eigen values and if this does not work, it uses a classical dichotomy
% technique.

Umin = evalseries(U0,0,sys.order);
Jacobian_str = Stab_Jacobian(sys,Umin);
[eig_val,eig_vec] = sys.Stability(Umin,Jacobian_str);
Eigen_init.values = eig_val;
Eigen_init.vectors = eig_vec;

nb_realpos1 = sum(real(Eigen_init.values)>sys.StabTol);
status_init = (nb_realpos1==0);
if status_init == 1
    status{1} = 'stable';
else
    status{1} = 'unstable';
end

Umax = evalseries(U0,Amax,sys.order);
Jacobian_str = Stab_Jacobian(sys,Umax);
[eig_val,eig_vec] = sys.Stability(Umax,Jacobian_str);
Eigen_end.values = eig_val;
Eigen_end.vectors = eig_vec;
nb_realpos2 = sum(real(eig_val)>sys.StabTol);
flag2 = (nb_realpos2==0);


Umid = Umax;
Amid = Amax;

    
if ((nb_realpos1 ~= nb_realpos2) && (sum(sum(isnan(eig_vec)) == 0)))
     disp('Change of the stability detected. Computation of the exact point.');
%     A1 = 0;
%     A2 = Amax;
%     nb_iter = 0;
%     
%     
%    % Tolerance on relative position of Achange :
%    A_tol = 1e-3;
%
%     try
%         warning('off','MATLAB:polyfit:RepeatedPointsOrRescale');
%         %% Polynomial fit - faster but sometimes does not work.
%         Alist = [0 Amax];
%         [~,min_ind1] = min(abs(real(Eigen_init.values)));
%         [~,min_ind2] = min(abs(real(Eigen_end.values)));
%         Vlist = [real(Eigen_init.values(min_ind1)) real(Eigen_end.values(min_ind2))];
%         
%         while abs(Alist(end)-Alist(end-1))/Amax>A_tol && nb_realpos1~=nb_realpos2 && nb_iter < sys.NRitemax
%             
%             Polynom = polyfit(Alist,Vlist,nb_iter+1);
%             Roots_P = roots(Polynom);
%             Roots_P = Roots_P(imag(Roots_P)==0);
%             Amid = Roots_P((Roots_P<=Amax & Roots_P>=0));
%             if numel(Amid) > 1
%                 error('Polynomial interpolation does not work in this case. Dichotomy used instead.');
%             end
%             Umid = evalseries(U0,Amid,sys.order);
%             Jacobian_str=Stab_Jacobian(sys,Umid);
%             
%             [eig_val,eig_vec] = sys.Stability(Umid,Jacobian_str);
%             nb_realpos = sum(real(eig_val)>sys.StabTol);
%             
%             if (nb_realpos == nb_realpos1) || (nb_realpos == nb_realpos2)
%                 Alist = [Alist Amid];
%                 [~,min_ind] = min(abs(real(eig_val)));
%                 Vlist = [Vlist real(eig_val(min_ind))];
%             else
%                 nb_iter = sys.NRitemax;
%             end
%             
%             if (Amid == 0) || (Amid == Amax)
%                 error('Polynomial interpolation does not work in this case. Dichotomy used instead.');
%             end
%             nb_iter = nb_iter+1;
%         end
%         
%         %disp(['Polynomial fit : nb_iter = ' num2str(nb_iter)]);
%         
%     catch
        
        A1 = 0;
        A2 = Amax;
        nb_iter = -1;
        
        %% Dichotomy - slower but always works
        while (nb_realpos1~=nb_realpos2 && nb_iter <= sys.NRitemax)
            
            Amid = (A2+A1)/2;
            Umid = evalseries(U0,Amid,sys.order);
            Jacobian_str=Stab_Jacobian(sys,Umid);
            
            [eig_val,eig_vec] = sys.Stability(Umid,Jacobian_str);
            nb_realpos = sum(real(eig_val)>sys.StabTol);
            
            if nb_realpos == nb_realpos1
                A1 = Amid;
                nb_realpos1 = nb_realpos;
            elseif nb_realpos == nb_realpos2
                A2 = Amid;
                nb_realpos2 = nb_realpos;
            else
                disp('StabilityComputation : Dichotomy did not converge.');
                nb_iter = sys.NRitemax;
            end
            nb_iter = nb_iter+1;
        end
        
        %disp(['Dichotomy : nb_iter = ' num2str(nb_iter)]);
        
%    end
    
    if nb_iter == sys.NRitemax
        msgbox({'@Syst/StabilityComputation : ' ; 'The location of the bifurcation is probably not accurate.'; 'Consider reducing the maximum value of Amax to compute it more accurately.'})
    end
    
    %% Post-processing
    if flag2 == 1
        status{2} = 'stable';
    else
        status{2} = 'unstable';
    end
    
    Uchange = Umid;
    Achange = Amid;
    
    bif_type = 'nothing';
    
    % Eigenvalues/Floquet exponents with positive real part
    eig1_pos = Eigen_init.values(real(Eigen_init.values)>sys.StabTol);
    eig2_pos = Eigen_end.values(real(Eigen_end.values)>sys.StabTol);
    nb_change = abs(numel(eig1_pos)-numel(eig2_pos));
    
    switch sys.type
        case {'HBQ','Q','HBM'}
            
            if nb_change == 2
                disp('Neimark-Sacker bifurcation detected.');
                bif_type = 'NS';
            elseif nb_change == 1
                % Floquet multipliers at the bifurcation point.
                Mult = exp(2*pi/Uchange(sys.neq)*eig_val);
                Mult_m1 = min(abs(Mult-1));
                Mult_p1 = min(abs(Mult+1));
                if Mult_m1 < Mult_p1
                    disp('Simple bifurcation detected.');
                    bif_type = 'B';
                else
                    disp('Period doubling bifurcation detected.');
                    bif_type = 'PD';
                end
            else
                disp('StabilityComputation : There could be more than one bifurcation. Consider reducing Amax.');
                bif_type = 'nothing';
            end
        case {'AQ','Equ'}
            
            if nb_change == 2
                switch sys.writing
                    case 'Hopf'
                        disp('Double-Hopf bifurcation detected.');
                        bif_type = 'HH';
                    case 'B'
                        disp('Zero-Hopf bifurcation detected.');
                        bif_type = 'ZH';
                    otherwise
                        disp('Hopf bifurcation detected.');
                        bif_type = 'H';
                end
            elseif nb_change == 1
                switch sys.writing
                    case 'Hopf'
                        disp('Zero-Hopf bifurcation detected.');
                        bif_type = 'ZH';
                    case 'B'
                        disp('Bogdanov-Takens bifurcation or Branching Point detected.');
                        bif_type = 'BT';
                    otherwise
                        disp('Simple bifurcation detected.');
                        bif_type = 'B';
                end
            else
                disp('StabilityComputation : There could be more than one bifurcation. Consider reducing Amax.');
                bif_type = 'nothing';
            end
        otherwise
            disp('It seems that the Stability Computation of this kind of system has not been implemented yet. The results could be wrong.');
    end
    
    Eigen.type = bif_type;
    
else
    Eigen.type = 'nothing';
    Uchange = 0;
    Achange = 0;
    if status_init == 1
        status{2} = 'stable';
    else
        status{2} = 'unstable';
    end
end

[~,ind] = sort(abs(real(eig_val)));
Eigen.values = eig_val(ind);
Eigen.vectors = eig_vec(:,ind);

StabData.Eigen_init = Eigen_init;
StabData.Eigen_end = Eigen_end;
StabData.Uchange = Uchange;
StabData.Achange = Achange;
StabData.Eigen = Eigen;
StabData.status = status;

end

function [Jacobian_str] = Stab_Jacobian(obj,U)
Jacobian_str = obj.Jacobian(U);
[Jacobian_str.LK_aux,Jacobian_str.UK_aux,Jacobian_str.pK_aux,Jacobian_str.qK_aux] = lu(Jacobian_str.dRauxdUaux,'vector');
Jacobian_str.Kfull = Jacobian_str.dRdU - ...
    Jacobian_str.dRdUaux(:,Jacobian_str.qK_aux)*(Jacobian_str.UK_aux\(Jacobian_str.LK_aux\Jacobian_str.dRauxdU(Jacobian_str.pK_aux,:)));
Jacobian_str.K = Jacobian_str.Kfull(:,1:end-1);
end





%%% Using test functions only.
%     switch obj.type
%         case 'HBQ'
%             floq_mult = exp(2*pi/Umid(obj.neq)*eig_val); % exp(T*rho)
%
%             prod_m1 = abs(prod(floq_mult-1));
%             prod_p1 = abs(prod(floq_mult+1));
%             for i=1:numel(floq_mult)-1
%                 bialt_prod(i) = prod( floq_mult(i)*floq_mult(i+1:end)-1);
%             end
%             min_abs_bialt_prod = min(abs(bialt_prod));
%
%             if prod_m1 < min(prod_p1,min_abs_bialt_prod)
%                 if abs(prod(floq_mult-1)) < A_tol
%                     disp('Simple bifurcation detected.');
%                     bif_type = 'SB';
%                 end
%             elseif prod_p1 < min_abs_bialt_prod
%                 if abs(prod(floq_mult+1)) < A_tol
%                     disp('Period doubling bifurcation detected.');
%                     bif_type = 'PD';
%                 end
%             else
%                 if min_abs_bialt_prod < A_tol
%                     disp('Neimark-Sacker bifurcation detected.');
%                     bif_type = 'NS';
%                 end
%             end
%         case 'AQ'
%
%             product_eig = abs(prod(eig_val));
%             prod_sum = Inf;
%             for i=1:numel(eig_val)-1
%                 prod_sum(i) = prod( eig_val(i) + eig_val(i+1:end));
%             end
%             min_abs_prodsum = min(abs(prod_sum));
%
%             if product_eig < min_abs_prodsum
%                 if product_eig < A_tol
%                     disp('Simple bifurcation detected');
%                     bif_type = 'SB';
%                 end
%             else
%                 if min_abs_prodsum < A_tol
%                     disp('Hopf bifurcation detected.');
%                     bif_type = 'Hopf';
%                 end
%             end
%     end
