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


function Hand = plotbranchHBM(sys,Sec_Diag,Idisp,Hdisp,bifpara_str)

Nb=length(Idisp);

% Color for successive plots.
Col = dispcolorML;

switch bifpara_str
    case 'omega'
        ind_bifpara = sys.neq;
    case 'lambda'
        ind_bifpara = sys.neq+1;
end


if isa(Sec_Diag,'CheckPoint')
    [Hand,Handleg,Leg] = plotsectionHBM(sys,Sec_Diag,Idisp,Hdisp);
else
    cellfun(@(sec)plotsectionHBM(sys,sec,Idisp,Hdisp),Sec_Diag(1:end-1),'uniformoutput',0);
    [Hand,Handleg,Leg] = plotsectionHBM(sys,Sec_Diag{end},Idisp,Hdisp);
end


legend(Handleg,Leg,'location','best'); % Slows the plotting procedure

ylabel('Harmonics amplitude')
xlabel(['$\' bifpara_str '$'])



    function [Hand,Handleg,Leg] = plotsectionHBM(sys,Section,Idisp,Hdisp)
        
        Hand=[];
        Handleg=[];
        
        st1 = Section.drawtype_init;
        st2 = Section.drawtype_end;
        
        bifpara = Section.Upp(ind_bifpara,:);
        
        H = sys.H;
        DHp1 = 2*H+1;
        
        
        %%% Test if there is a bifurcation in the Section
        if ~strcmp(Section.Eigen.type,'nothing')
            ind_change = Section.ind_change-1;
            
            bifparastab = Section.Ustab(ind_bifpara);
            
            bifpara = [bifpara(1:ind_change-1),bifparastab,bifpara(ind_change:end)];
        end
        
        for ii=1:Nb
            idisp = Idisp(ii);
            hdisp = Hdisp(ii);
            
            cdisp = Col(ii,:);
            
            if idisp <= sys.nz
                I0   = (idisp-1)*DHp1+1;
                Icos = (idisp-1)*DHp1+1+hdisp;
                Isin = (idisp-1)*DHp1+1+H+hdisp;
            else
                I0   = sys.neq+1+ (idisp-1-sys.nz)*DHp1+1;
                Icos = sys.neq+1+ (idisp-1-sys.nz)*DHp1+1+hdisp;
                Isin = sys.neq+1+ (idisp-1-sys.nz)*DHp1+1+H+hdisp;
            end
            
            Hampl   = sqrt(Section.Upp(Icos,:).^2+Section.Upp(Isin,:).^2);
            
            %%% Test if there is a bifurcation in the Section
            if strcmp(Section.Eigen.type,'nothing')
                hand = plot(bifpara,Hampl,st1,'Color',cdisp,'linewidth',1);
                Hand = [Hand;hand];
            else
                Hampl_stab   = sqrt(Section.Ustab(Icos,:).^2+Section.Ustab(Isin,:).^2);
                Hampl = [Hampl(1:ind_change-1),Hampl_stab,Hampl(ind_change:end)];
                
                hand = plot(bifpara(1:ind_change),Hampl(1:ind_change),st1,'Color',cdisp,'linewidth',1);
                Hand = [Hand;hand];
                hand = plot(bifpara(ind_change:end),Hampl(ind_change:end),st2,'Color',cdisp,'linewidth',1);
                Hand = [Hand;hand];
                hand = plot(bifpara(ind_change),Hampl(ind_change),'.','Color',cdisp,'linewidth',1);
                Hand = [Hand;hand];
                
            end
            
            Handleg(ii) = plot(NaN,NaN,'-','Color',cdisp,'linewidth',1);
            Leg{ii}=['z$_{',num2str(idisp),'}$H',num2str(hdisp)];
        end
        
    end


end