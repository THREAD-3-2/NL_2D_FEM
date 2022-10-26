function [EigenVal,FloqMult] = calcdiagStab(sys,Diagram)
% function [EigenVal] = calcdiagStab(sys,Diagram)
%
% From the Diagram structure, computes all the Eigen information at all the
% points of discretization of the series and puts them in EigenVal matrix.
%

nb_Section = length(Diagram);
nb_pts = size(Diagram{1}.Upp,2);
nb_pts_diag = nb_Section*nb_pts;

EigenVal = zeros(sys.nz,nb_pts_diag);
FloqMult = zeros(size(EigenVal));

index = 1;

UUend=Diagram{1}.Upp(:,1);

for isec = 1:nb_Section
    for ipts = 1:nb_pts
        Uii = Diagram{isec}.Upp(:,ipts);
        
        eigval = get_eigval(sys,Uii);
        
        %%% Bifurcation ...
        if numel(Diagram{isec}.Ustab) > 1
            eigval0 = EigenVal(:,index-1);
            nb_pos0 = sum(real(eigval0)>sys.StabTol);
            nb_pos = sum(real(eigval)>sys.StabTol);
            diff_pos = abs(nb_pos0-nb_pos);
            % and eigenvalues at bifurcation.
            if diff_pos > 0
                eigvalBIF = Diagram{isec}.Eigen.values; 
                
                %%% Put the eigenvalue in the full matrix.
                nb_eigval = numel(eigvalBIF);
                EigenVal(1:nb_eigval,index) = eigvalBIF;
                FloqMult(1:nb_eigval,index) = exp(2*pi/Uii(sys.neq)*eigvalBIF);
                if nb_eigval < sys.nz
                    EigenVal(nb_eigval+1:end,index) = NaN;
                    FloqMult(nb_eigval+1:end,index) = NaN;
                end
                index = index+1;
            end
        end
        
        %%% Put the eigenvalue in the full matrix.
        nb_eigval = numel(eigval);
        EigenVal(1:nb_eigval,index) = eigval;
        FloqMult(1:nb_eigval,index) = exp(2*pi/Uii(sys.neq)*eigval);
        if nb_eigval < sys.nz
            EigenVal(nb_eigval+1:end,index) = NaN;
            FloqMult(nb_eigval+1:end,index) = NaN;
        end
    
        index = index+1;
    
    end
    
    UU0= Diagram{isec}.Upp(:,1);
    testjump= norm(UUend-UU0)>0; % test for a jump in the diagram
    if testjump
        EigenVal(:,index+1) = NaN;
        FloqMult(:,index+1) = NaN;
        index = index+1;
    end
    UUend= Diagram{isec}.Upp(:,end);

end


end


function eigval = get_eigval(sys,U)
%%% Computation of the Jacobian matrix at point Uii
Jacobian = sys.Jacobian(U);
%%% Condensation of the auxiliary variables
[Jacobian.LK_aux,Jacobian.UK_aux,Jacobian.pK_aux,Jacobian.qK_aux] = lu(Jacobian.dRauxdUaux,'vector');
Jacobian.Kfull = Jacobian.dRdU - ...
    Jacobian.dRdUaux(:,Jacobian.qK_aux)*(Jacobian.UK_aux\(Jacobian.LK_aux\Jacobian.dRauxdU(Jacobian.pK_aux,:)));
%%% Remove the column dR/dlambda
Jacobian.K = Jacobian.Kfull(:,1:end-1);

%%% Computation of the eigenvalues
eigval = sys.Stability(U,Jacobian);
eigval = sort(eigval);
end