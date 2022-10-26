% function Hand = plotdiagnostabHBM(sys,Diag,Idisp,Hdisp,bifpara_str)
%    plots the amplitude of the harmonics Hdisp of the variables Idisp computed
%    by the Manlab/HBM method as a function of the bifurcation parameter
%    specified by bifpara_str ( = 'omega' or 'lambda').
%
%    The input Diag can be:
%    - a CheckPoint object 
%        (i.e. obtained by the Export Section button of the Manlab interface) 
%    - a cell array of CheckPoint 
%        (i.e. obtained by the Export Diagram button of the Manlab interface)
%    - or a Diag structure more convenient for ploting bifurcation diagrams
%        (obtained by calcdiagUpp.m function)
%        
%  
%
%    By L. Guillot and O. Thomas / Nov. 2018 - Mars 2019


function Hand = plotdiagnostabHBM(sys,Diag,Idisp,Hdisp,bifpara_str)

if ~isstruct(Diag)
    [Diag] = calcdiagUpp(sys,Diag);
end

DiagUpp = Diag.DiagUpp;

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
    hdisp = Hdisp(ii);
    cdisp = Col(ii,:);
    
    if idisp <= sys.nz
        if hdisp == 0
            I0 = (idisp-1)*DHp1+1;
            Norm = sqrt(DiagUpp(I0,:).^2);
        else
            I0 = 0;
            Icos = (idisp-1)*DHp1+1+hdisp;
            Isin = (idisp-1)*DHp1+1+H+hdisp;
            Norm = sqrt(DiagUpp(Icos,:).^2+DiagUpp(Isin,:).^2);
        end
    else
        if hdisp == 0
            I0 = sys.neq+1+ (idisp-1-sys.nz)*DHp1+1;
            Norm = sqrt(DiagUpp(I0,:).^2);
        else
            I0 = 0;
            Icos = sys.neq+1+ (idisp-1-sys.nz)*DHp1+1+hdisp;
            Isin = sys.neq+1+ (idisp-1-sys.nz)*DHp1+1+H+hdisp;
            Norm = sqrt(DiagUpp(Icos,:).^2+DiagUpp(Isin,:).^2);
        end
    end
        
    hand = plot(bifpara,Norm,'Color',cdisp,'linewidth',1);
    Hand = [Hand;hand];
    Handleg(ii) = plot(NaN,NaN,'-','Color',cdisp);
    Leg{ii}=['z_',num2str(idisp),' H_',num2str(hdisp)];
end

legend(Handleg,Leg,'location','best'); % Slows the plotting procedure

ylabel('harmonic amplitude')
xlabel(['\' bifpara_str])


end