% plotbranchampphaseHBM   Branch plot for Manlab/HBM.
%    Hand = plotbranchampphaseHBM(sys,Sec_Diag,Idisp,Hdisp,bifpara_str) plots the
%    amplitude and the phase difference of the harmonics components
%    of the solution computed by the Manlab/HBM method
%    as a function of the angular frequency (omega=U(end,:)).
%
%    sys is an object containing the informations about the system solved.
%
%    Sec_Diag is an object containing a part of a bifurcation diagram
%    computed with Manlab.
%
%    The displayed amplitude zi and phase phii of the i-th harmonic
%    component are defined by:
%
%    z(t) = z0 + z1 cos(omega t + phi1) + ... + zi cos(i*omega t + phii) + ...
%
%    Idisp contains the indices of the entries of the initial physical
%    unknown vector u(t) to be displayed and
%    Hdisp the corresponding harmonics number.
%
%    fig is a figure handle to plot the bargraph on.
%
%    It returns in Hand a column vector of handles to lineseries objects,
%    one handle per plotted line.
%
%    Example: Hand=plotbranchampphaseHBM(sys,Section,[1 1 1],[1 3 5],1)
%             plots the three first odd harmonics (H=1, 3, 5) of u1.
%
%    By L. Guillot and O. Thomas / Nov. 2018

function Hand = plotbranchampphaseHBM(sys,Sec_Diag,Idisp,Hdisp,bifpara_str)

subplot(2,1,1)
hold on
subplot(2,1,2)
hold on

Nb=length(Idisp);

% Color for successive plots.
Col = dispcolorML;

switch bifpara_str
    case 'omega'
        ind_bifpara = sys.neq;
    case 'lambda'
        ind_bifpara = sys.neq+1;
end

H = sys.H;
DHp1 = 2*H+1;

if isa(Sec_Diag,'CheckPoint')
    [Hand,Handleg,Leg] = plotsectionHBM(sys,Sec_Diag,Idisp,Hdisp);
else
    cellfun(@(sec)plotsectionHBM(sys,sec,Idisp,Hdisp),Sec_Diag(1:end-1),'uniformoutput',0);
    [Hand,Handleg,Leg] = plotsectionHBM(sys,Sec_Diag{end},Idisp,Hdisp);
end

subplot(2,1,1)
legend(Handleg,Leg,'location','best') % Slows the plotting procedure
box on
ylabel('Harmonics amplitude')
xlabel(['\' bifpara_str])
subplot(2,1,2)
box on
ylabel('Harmonics phase [\pi rad]')
xlabel(['\' bifpara_str])



    function [Hand,Handleg,Leg] = plotsectionHBM(sys,Section,Idisp,Hdisp)
        
        st1 = Section.drawtype_init;
        st2 = Section.drawtype_end;
        
        Hand=[];
        Handleg=[];
        
        bifpara = Section.Upp(ind_bifpara,:);
        
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
            
            if hdisp == 0
                Amp   = Section.Upp(I0,:);
                Phase = NaN*ones(size(Amp));
            else
                Amp   = sqrt(Section.Upp(Icos,:).^2+Section.Upp(Isin,:).^2);
                Phase = -atan2(Section.Upp(Isin,:),Section.Upp(Icos,:)) / pi;
            end
            
            if strcmp(Section.Eigen.type,'nothing')
                subplot(2,1,1)
                hand = plot(bifpara,Amp,st1,'Color',cdisp,'linewidth',1);
                Hand = [Hand;hand];
                subplot(2,1,2)
                hand = plot(bifpara,Phase,st1,'Color',cdisp,'linewidth',1);
                Hand = [Hand;hand];
            else
                if hdisp == 0
                    Amp   = Section.Ustab(I0,:);
                    Phase = NaN*ones(size(Amp));
                else
                    Ampstab   = sqrt(Section.Ustab(Icos,:).^2+Section.Ustab(Isin,:).^2);
                    Phasestab = -atan2(Section.Ustab(Isin,:),Section.Ustab(Icos,:)) / pi;
                    
                    Amp = [Amp(1:ind_change-1),Ampstab,Amp(ind_change:end)];
                    Phase = [Phase(1:ind_change-1),Phasestab,Phase(ind_change:end)];
                end
                
                subplot(2,1,1)
                hand = plot(bifpara(1:ind_change),Amp(1:ind_change),st1,'Color',cdisp,'linewidth',1);
                Hand = [Hand;hand];
                hand = plot(bifpara(ind_change:end),Amp(ind_change:end),st2,'Color',cdisp,'linewidth',1);
                Hand = [Hand;hand];
                hand = plot(bifpara(ind_change),Amp(ind_change),'.','Color',cdisp,'linewidth',1);
                Hand = [Hand;hand];
                
                subplot(2,1,2)
                hand = plot(bifpara(1:ind_change),Phase(1:ind_change),st1,'Color',cdisp,'linewidth',1);
                Hand = [Hand;hand];
                hand = plot(bifpara(ind_change:end),Phase(ind_change:end),st2,'Color',cdisp,'linewidth',1);
                Hand = [Hand;hand];
                hand = plot(bifpara(ind_change),Phase(ind_change),'.','Color',cdisp,'linewidth',1);
                Hand = [Hand;hand];
                
            end
            Handleg(ii) = plot(NaN,NaN,'-','Color',cdisp,'linewidth',1);
            Leg{ii}=['H_',num2str(hdisp),'z_',num2str(idisp)];
        end
        
    end


end

