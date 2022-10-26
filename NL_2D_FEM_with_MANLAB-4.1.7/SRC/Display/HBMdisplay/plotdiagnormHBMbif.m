% function Hand = plotbranchnormHBM(sys,Diag,Idisp,bifpara_str)
%    plots the L2 norm of the variables Idisp computed
%    by the Manlab/HBM method as a function of the bifurcation parameter.
%
%    It also displays the bifurcation points:
%    - B for a simple bifurcation (saddle-node, pitchfork...)
%    - PD for period doubling bifurcation
%    - NS for Neimark-Sacker bifurcation
%
%    By L. Guillot and O. Thomas / Nov. 2018

function Hand = plotdiagnormHBMbif(sys,Diag,Idisp,bifpara_str)

if ~isstruct(Diag)
    [Diag] = calcdiagUpp(sys,Diag);
end

DiagUpp = Diag.DiagUpp;
Stabinfo = Diag.Stabinfo;
change = Diag.change;

Nb=length(Idisp);
Hand=[];
Handleg=[];

% Color for successive plots.
Col = dispcolorML;

switch bifpara_str
    case 'omega'
        ind_bifpara = sys.neq;
    case 'lambda'
        ind_bifpara = sys.neq+1;
end

bifpara = DiagUpp(ind_bifpara,:);

strstab = '-';
strunstab = ':';

H = sys.H;
DHp1 = 2*H+1;

for ii=1:Nb
    idisp = Idisp(ii);
    
    cdisp = Col(ii,:);
    
    if idisp <= sys.nz
        I0   = (idisp-1)*DHp1+1;
        Icos = (idisp-1)*DHp1+1+(1:H);
        Isin = (idisp-1)*DHp1+1+H+(1:H);
    else
        I0   = sys.neq+1+ (idisp-1-sys.nz)*DHp1+1;
        Icos = sys.neq+1+ (idisp-1-sys.nz)*DHp1+1+(1:H);
        Isin = sys.neq+1+ (idisp-1-sys.nz)*DHp1+1+H+(1:H);
    end
    
    L2norm   = sqrt( DiagUpp(I0,:).^2 + .5*sum(DiagUpp(Icos,:).^2+DiagUpp(Isin,:).^2,1));
    
    if numel(change) == 0
        
        if Stabinfo(1) == 1; strplot = strstab; else; strplot = strunstab; end
        hand = plot(bifpara,L2norm,strplot,'Color',cdisp);
        Hand = [Hand;hand];
        
    else
        
        for istab=1:numel(change)
            ind_change = change(istab);
            if istab == 1; indplot = 1:ind_change; else; indplot = change(istab-1):ind_change; end
            if Stabinfo(indplot(2)) == 1; strplot = strstab; else; strplot = strunstab; end
            
            L2normstab   = sqrt( DiagUpp(I0,ind_change).^2 + .5*sum(DiagUpp(Icos,ind_change).^2+DiagUpp(Isin,ind_change).^2,1));
            bifparastab = DiagUpp(ind_bifpara,ind_change);
            
            hand = plot(bifpara(indplot),L2norm(indplot),strplot,'Color',cdisp);
            Hand = [Hand;hand];
            hand = plot(bifpara(ind_change),L2norm(ind_change),'p','Color',cdisp);
            Hand = [Hand;hand];
            
            text(bifparastab,L2normstab,char(Stabinfo(ind_change)),'margin',3,'verticalalignment','Bottom','fontsize',18)
            box on;
        end
        
        if Stabinfo(ind_change+1) == 1; strplot = strstab; else; strplot = strunstab; end
        hand = plot(bifpara(ind_change:end),L2norm(ind_change:end),strplot,'Color',cdisp);
        Hand = [Hand;hand];
        
    end
    
    Handleg(ii) = plot(NaN,NaN,'-','Color',cdisp);
    Leg{ii}=['z_',num2str(idisp)];
end

legend(Handleg,Leg,'location','best'); % Slows the plotting procedure

ylabel('L^2-norm')
xlabel(['\' bifpara_str])


end