function Hand = plotperiodHBM(sys,U,Idisp,time)
%    Hand = plotperiodHBM(sys,U,Idisp,time) 
%    
%    Plots the periodic solution in the time domain, from the Fourier
%    expansion of the variables Icalc. 
%
%    U is the vector of unknowns of the final algebraic system solved
%    by Manlab.
%
%    Idisp contains the indices of the entries of the initial physical
%    unknown vector u(t) to be displayed. It must be of length 2.
%
%    It returns in Hand a column vector of handles to lineseries objects,
%    one handle per plotted line.
%
%    By L. Guillot and O. Thomas / Jan. 2022


Hand = [];

H = sys.H;

if nargin < 4
  time=linspace(0,1,100*H);
end

Utime = calcperiodHBM(sys,U,Idisp,time);

Col = dispcolorML;
Leg = cell(length(Idisp),1);

for ii=1:length(Idisp)
    idisp = Idisp(ii);
    cdisp = Col(ii,:);
    hand=plot(time,Utime(:,ii),'Color',cdisp);
    Hand=[Hand;hand];
    Leg{ii} = ['z_{' num2str(idisp) '}'];
end
xlabel('dimensionless time');
legend(Hand,Leg);
box on

