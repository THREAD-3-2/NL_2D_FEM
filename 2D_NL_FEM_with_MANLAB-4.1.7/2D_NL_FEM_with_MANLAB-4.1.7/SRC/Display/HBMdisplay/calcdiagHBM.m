% calcbranchHBM   Diagram computation for Manlab/HBM.
% function Hampl = calcdiagHBM(sys,Diag,Icalc,Hcalc)
%    Hampl : amplitude of the harmonics components of the solution computed by the
%    Manlab/HBM.
%
%    sys is an object containing the informations about the system solved.
%
%    Diag is an object containing a part of a bifurcation diagram
%    computed with Manlab.
%
%    Icalc contains the indices of the entries of the initial physical
%    unknown vector z(t) to be computed and Hcalc the corresponding
%    harmonics number.
%
%    It returns Hampl, a matrix of size (length(Icalc),nb_pts),
%    containing the values of the branch amplitude where nb_pts is the
%    number of points to plot the projected bifurcation diagram.
%
%    By L. Guillot and O. Thomas / Nov. 2018

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
%    By L. Guillot and O. Thomas / Nov. 2018 - April 2019


function Hampl = calcdiagHBM(sys,Diag,Icalc,Hcalc)

if ~isstruct(Diag)
    [Diag] = calcdiagUpp(sys,Diag);
end

DiagUpp = Diag.DiagUpp;

Nb=length(Icalc);
Hand=[];
Handleg=[];

H = sys.H;
DHp1 = 2*H+1;

Hampl = zeros(Nb,size(DiagUpp,2));

for ii=1:Nb
    idisp = Icalc(ii);
    hdisp = Hcalc(ii);
    
    if idisp <= sys.nz
        if hdisp == 0
            I0 = (idisp-1)*DHp1+1;
            Hampl(ii,:) =  abs(DiagUpp(I0,:));
        else
            Icos = (idisp-1)*DHp1+1+hdisp;
            Isin = (idisp-1)*DHp1+1+H+hdisp;
            Hampl(ii,:) = sqrt(DiagUpp(Icos,:).^2+DiagUpp(Isin,:).^2);
        end
    else
        if hdisp == 0
            I0 = sys.neq+1+ (idisp-1-sys.nz)*DHp1+1;
            Hampl(ii,:) =  abs(DiagUpp(I0,:));
        else
            Icos = sys.neq+1+ (idisp-1-sys.nz)*DHp1+1+hdisp;
            Isin = sys.neq+1+ (idisp-1-sys.nz)*DHp1+1+H+hdisp;
            Hampl(ii,:) = sqrt(DiagUpp(Icos,:).^2+DiagUpp(Isin,:).^2);
        end
    end 
        
end

end