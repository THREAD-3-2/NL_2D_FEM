function [Utime] = calctimeQPHBM(sys,U,Icalc,time)
% function [Utime] = calctimeQPHBM(sys,U,Icalc,time)
%
% Compute the periodic solution in the time domain, from the Fourier
% developments of the variables Icalc.
%
% time is the grid on which the solution should be represented. If not
% specified, the solution is computed over one period.
% The time is dimensionless : 1 is the period.

[Z,omega1,omega2]=sys.get_Ztot(U);
if nargin < 4
    Tmoy = pi/omega1 + pi/omega2;
    time=linspace(0,20*Tmoy,1e4);
end
time = time(:);
nt=length(time);

H = sys.H;
%coef_per_var = (H+1)*4*H+1;

D1vec = [zeros(H,1) ; reshape((1:H).*ones(2*H+1,1),2*H*H+H,1)]';
D2vec = [(1:H)' ; reshape((-H:1:H)'.*ones(1,H),2*H*H+H,1)]';

VectCos=cos(2*pi*time*(D1vec*omega1 + D2vec*omega2));
VectSin=sin(2*pi*time*(D1vec*omega1 + D2vec*omega2));

MAT=[ ones(nt,1) , VectCos , VectSin ];  % MAT(nt,coef_per_var)
Utime=MAT*Z(:,Icalc);  % Ut(nt,nz)


end

