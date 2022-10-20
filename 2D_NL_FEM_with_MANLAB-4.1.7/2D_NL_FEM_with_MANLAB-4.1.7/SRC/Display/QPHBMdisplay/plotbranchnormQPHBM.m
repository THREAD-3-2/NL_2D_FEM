% function Hand = plotbranchnormHBM(sys,Sec_Diag,Idisp,bifpara_str)
%    plots the L2 norm of the variables Idisp computed
%    by the Manlab/HBM method as a function of the bifurcation parameter.
%
%    By L. Guillot - Mars 2019

function Hand = plotbranchnormQPHBM(sys,Sec_Diag,Idisp,bifpara_str)

Nb=length(Idisp);

% Color for successive plots.
Col = dispcolorML;

switch bifpara_str
    case 'omega1'
        ind_bifpara = sys.neq-1;
    case 'omega2'
        ind_bifpara = sys.neq;
    case 'lambda'
        ind_bifpara = sys.neq+1;
end

if isa(Sec_Diag,'CheckPoint')
    [Hand,Handleg,Leg] = plotsectionHBM(sys,Sec_Diag,Idisp);
else
    cellfun(@(sec)plotsectionHBM(sys,sec,Idisp),Sec_Diag(1:end-1),'uniformoutput',0);
    [Hand,Handleg,Leg] = plotsectionHBM(sys,Sec_Diag{end},Idisp);
end

legend(Handleg,Leg,'location','best'); % Slows the plotting procedure

ylabel('L^2-norm')
xlabel(['\' bifpara_str])


    function [Hand,Handleg,Leg] = plotsectionHBM(sys,Section,Idisp)
        
        Hand=[];
        Handleg=[];
        
        st1 = Section.drawtype_init;
        st2 = Section.drawtype_end;
        
        bifpara = Section.Upp(ind_bifpara,:);
        
            H = sys.H;
            coef_per_var = (2*H+1)^2;
            nb_cosin = (coef_per_var-1)/2;

        if ~strcmp(st1,st2)
            ind_change = Section.ind_change-1;
            
            bifparastab = Section.Ustab(ind_bifpara);
            
            bifpara = [bifpara(1:ind_change-1),bifparastab,bifpara(ind_change:end)];
        end
        
        for ii=1:Nb
            idisp = Idisp(ii);
            
            cdisp = Col(ii,:);
            
            if idisp <= sys.nz
                I0   = (idisp-1)*coef_per_var+1;
                Icos = (idisp-1)*coef_per_var+1+(1:nb_cosin);
                Isin = (idisp-1)*coef_per_var+1+nb_cosin+(1:nb_cosin);
            else
                I0   = sys.neq+1+ (idisp-1-sys.nz)*coef_per_var+1;
                Icos = sys.neq+1+ (idisp-1-sys.nz)*coef_per_var+1+(1:nb_cosin);
                Isin = sys.neq+1+ (idisp-1-sys.nz)*coef_per_var+1+nb_cosin+(1:nb_cosin);
            end
            
            L2norm   = sqrt( Section.Upp(I0,:).^2 + .5*sum(Section.Upp(Icos,:).^2+Section.Upp(Isin,:).^2,1));
            
            if strcmp(st1,st2)
                hand = plot(bifpara,L2norm,st1,'Color',cdisp);
                Hand = [Hand;hand];
            else
                L2normstab   = sqrt( Section.Ustab(I0,:).^2 + .5*sum(Section.Ustab(Icos,:).^2+Section.Ustab(Isin,:).^2,1));
                L2norm = [L2norm(1:ind_change-1),L2normstab,L2norm(ind_change:end)];
                
                
                hand = plot(bifpara(1:ind_change),L2norm(1:ind_change),st1,'Color',cdisp);
                Hand = [Hand;hand];
                hand = plot(bifpara(ind_change:end),L2norm(ind_change:end),st2,'Color',cdisp);
                Hand = [Hand;hand];
                hand = plot(bifpara(ind_change),L2norm(ind_change),'.','Color',cdisp);
                Hand = [Hand;hand];
                
            end
            
            Handleg(ii) = plot(NaN,NaN,'-','Color',cdisp);
            Leg{ii}=['z_',num2str(idisp)];
        end
        
    end



end