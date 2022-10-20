function [Utime] = calcperiodHBM(sys,U,Icalc,time)
% function [Utime] = calcperiodHBM(sys,U,Icalc,time)
%
% Compute the periodic solution in the time domain, from the Fourier
% developments of the variables Icalc.
%
% time is the grid on which the solution should be represented. If not
% specified, the solution is computed over one period. 
% The time is dimensionless : 1 is the period.


switch sys.type
    case {'HBQ','HBM'}
        H = sys.H;
        Z=sys.get_Ztot(U);
        if nargin < 4
            time=linspace(0,1,100*H);
        end
        time = time(:);
        nt=length(time);
        
        VectCos=cos(2*pi*time*(1:H));
        VectSin=sin(2*pi*time*(1:H));
        
        MAT=[ ones(nt,1) , VectCos , VectSin ];  % MAT(nt,2*H+1)
        Utime=MAT*Z(:,Icalc);  % Ut(nt,nz)
    case 'HBQdiff'
        H = sys.H;
        Z=sys.get_Ztot(U);
        if nargin < 4
            time=linspace(0,1,100*H);
        end
        time = time(:);
        nt=length(time);
        Utime = zeros(nt,length(Icalc));
        for i=1:length(Icalc)
            ivar = Icalc(i);
            Hvar = H(ivar);
            VectCos=cos(2*pi*time*(1:Hvar));
            VectSin=sin(2*pi*time*(1:Hvar));
            
            MAT=[ ones(nt,1) , VectCos , VectSin ];  % MAT(nt,2*H+1)
            Utime(:,i)=MAT*Z{ivar};  % Ut(nt,nz)
        end
        
    case {'HBQ_QP','QPHBM'}
        disp('calcperiodHBM : Quasi-periodic solutions not available yet.');
end


end

