% function Hand = plotbranchnormHBM(sys,Section,Idisp,bifpara_str)
%    plots the L2 norm of the variables Idisp computed
%    by the Manlab/HBM method as a function of the bifurcation parameter.
%
%    It also displays the bifurcation points:
%    - B for a simple bifurcation (saddle-node, pitchfork...)
%    - PD for period doubling bifurcation
%    - NS for Neimark-Sacker bifurcation
%
%    By L. Guillot and O. Thomas / Nov. 2018

function Hand = plotdiagwaterfallHBM(sys,Diag,Idisp,Ia,bifpara_str)

if ~isstruct(Diag)
    [Diag] = calcdiagUpp(sys,Diag);
end

DiagUpp = Diag.DiagUpp;
Stabinfo = Diag.Stabinfo;
change = Diag.change;

if numel(Idisp)~=2
    disp('Two variables exactly must be specified in Idisp (third argument).');
end

Ndisp=length(Idisp);
Na = length(Ia);
Hand=[];

if min(Ia)>size(DiagUpp,2)
    disp('min(Ia) should be <= size(Diag.DiagUpp,2)')
    Ia=1:size(U,2);
end

if max(Ia)>size(DiagUpp,2)
    disp('max(Ia) should be <= size(Diag.DiagUpp,2)')
    Ia=1:size(U,2);
end

% Color for successive plots.
Col = dispcolorML(3);

switch bifpara_str
    case 'omega'
        ind_bifpara = sys.neq;
    case 'lambda'
        ind_bifpara = sys.neq+1;
end

bifpara = DiagUpp(ind_bifpara,:);
bifparaplot = bifpara(Ia);
stabplot = Stabinfo(Ia);

% Calcul des orbites périodiques
% ==============================
nech = 1000;
time = linspace(0,1,nech);

for jj=1:Na
    ia=Ia(jj);
    for ii=1:Ndisp
        idisp=Idisp(ii);
        utime(:,ii)=calcperiodHBM(sys,DiagUpp(:,ia),idisp,time);
    end
    
    if stabplot(jj) == 1; cplot = Col(1,:); linwid = 1.5;
    elseif stabplot(jj) == 0; cplot = Col(2,:); linwid = 1.5;
    else; cplot = Col(3,:); linwid = 3;
    end
    
    hand = plot3(bifparaplot(jj)*ones(nech,1),utime(:,1),utime(:,2),'color',cplot,'Linewidth',linwid);
    Hand = [Hand hand];
end

xlabel(['\' bifpara_str])
ylabel(['z_',num2str(Idisp(1))])
zlabel(['z_',num2str(Idisp(2))])
view(20,5)
box on

end