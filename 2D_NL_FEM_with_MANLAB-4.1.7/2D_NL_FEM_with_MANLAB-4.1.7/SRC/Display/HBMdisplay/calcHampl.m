function [Hampl] = calcHampl(sys,U,Icalc,Hcalc)
% function [Hampl] = calcHampl(sys,U,Icalc,Hcalc)
%
% For the variable "Icalc" of the original system,
% compute the amplitude of the harmonic "Hcalc".
%
% If "Hcalc" is not given, Hampl(:,Icalc(i)) gives all the harmonics
% amplitude of the variable Icalc(i).

Nb = numel(Icalc);

switch sys.type
    
    case {'HBQ','HBM'}
        [Z]=sys.get_Ztot(U);
        Zcompl = sys.real_to_compl(Z);
        
        if nargin < 4
            Hampl = zeros(sys.H+1,Nb);
                Hampl(:,Icalc) = 2*abs(Zcompl(:,Icalc));
                Hampl(1,:) = Hampl(1,:)/2;
        else
            Hampl = zeros(Nb,1);
            for i=1:Nb
                if Hcalc(i) == 0
                    factor = 1;
                else
                    factor = 2;
                end
                Hampl(i) = factor*abs(Zcompl(Hcalc(i)+1,Icalc(i)));
            end
        end
        
        
    case {'HBQ_QP','QPHBM'}
        disp('caclHampl : Quasi-periodic solutions not available yet.');
end


end

