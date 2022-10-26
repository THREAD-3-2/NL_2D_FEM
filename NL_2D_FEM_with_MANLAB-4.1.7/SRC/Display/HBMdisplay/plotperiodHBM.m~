function Hand = plotperiodHBM(sys,U,Idisp)
%    Hand = plotphasediagHBM(sys,U,Idisp) plots a phase diagram
%    of the periodic solution computed by the Manlab/HBM method.
%
%    U is the vector of unknowns of the final algebraic system solved
%    by Manlab.
%
%    Idisp contains the indices of the entries of the initial physical
%    unknown vector u(t) to be displayed. It must be of length 2.
%
%    fig is a figure/axis handle to plot the signal on.
%
%    It returns in Hand a column vector of handles to lineseries objects,
%    one handle per plotted line.
%
%    Example: Hand=plotphasediagHBM(sys,U,[1 2],1)
%
%    By L. Guillot and O. Thomas / Nov. 2018


Hand = [];

H = sys.H;

time = linspace(0,1,100*H);
Utime = calcperiodHBM(sys,U,Idisp);

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

end

