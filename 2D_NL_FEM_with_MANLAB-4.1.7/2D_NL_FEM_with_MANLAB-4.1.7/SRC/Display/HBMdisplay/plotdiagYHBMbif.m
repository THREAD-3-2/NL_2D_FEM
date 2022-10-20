% function Hand = plotbranchHBM(sys,Section,Idisp,Hdisp,bifpara_str)
%    plots the amplitude of the harmonics Hdisp of the variables Idisp computed
%    by the Manlab/HBM method as a function of the bifurcation parameter
%    specified by bifpara_str ( = 'omega' or 'lambda').
%
%    It also displays the bifurcation points:
%    - B for a simple bifurcation (saddle-node, pitchfork...)
%    - PD for period doubling bifurcation
%    - NS for Neimark-Sacker bifurcation
%
%
%    By L. Guillot and O. Thomas / Nov. 2018 - Mars 2019


function Hand = plotdiagYHBMbif(sys,Diag,Y,bifpara_str)

if ~isstruct(Diag)
    [Diag] = calcdiagUpp(sys,Diag);
end

DiagUpp = Diag.DiagUpp;
Stabinfo = Diag.Stabinfo;
change = Diag.change;

Nb=size(Y,1);
Hand=[];

% Color for successive plots.
Col = dispcolorML2(Nb);

switch bifpara_str
    case 'omega'
        ind_bifpara = sys.neq;
    case 'lambda'
        ind_bifpara = sys.neq+1;
end


bifpara = DiagUpp(ind_bifpara,:);

strstab = '-';
strunstab = ':';

for ii=1:Nb
    
    cdisp = Col(ii,:);
    Ydisp = Y(ii,:);
    if numel(change) == 0
        
        if Stabinfo(1) == 1; strplot = strstab; else; strplot = strunstab; end
        hand = plot(bifpara,Ydisp,strplot,'Color',cdisp);
        Hand = [Hand;hand];
        
    else
        
        if numel(change) == 0
            if Stabinfo(1) == 1; strplot = strstab; else; strplot = strunstab; end
            
            hand = plot(bifpara,Ydisp,strplot,'Color',cdisp);
            Hand = [Hand;hand];
            
        else
            for istab=1:numel(change)
                ind_change = change(istab);
                if istab == 1; indplot = 1:ind_change; else; indplot = change(istab-1):ind_change; end
                if Stabinfo(indplot(2)) == 1; strplot = strstab; else; strplot = strunstab; end
                
                Ystab = Ydisp(ind_change);
                
                bifparastab = DiagUpp(ind_bifpara,ind_change);
                
                hand = plot(bifpara(indplot),Ydisp(indplot),strplot,'Color',cdisp);
                Hand = [Hand;hand];
                hand = plot(bifparastab,Ydisp(ind_change),'.','Color',cdisp);
                Hand = [Hand;hand];
                
                text(bifparastab,Ystab,char(Stabinfo(ind_change)),'margin',3,'verticalalignment','Bottom')
                box on;
            end
            
            if Stabinfo(ind_change+1) == 1; strplot = strstab; else; strplot = strunstab; end
            hand = plot(bifpara(ind_change:end),Ydisp(ind_change:end),strplot,'Color',cdisp);
            Hand = [Hand;hand];
            
        end
    end
end


xlabel(['\' bifpara_str])


end