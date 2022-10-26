% calcbranchHBM   Branch computation for Manlab/HBM.
% function Hampl = calcbranchHBM(sys,Section,Icalc,Hcalc)
%    Hampl : amplitude of the harmonics components of the solution computed by the
%    Manlab/HBM.
%
%    sys is an object containing the informations about the system solved.
%
%    Section is an object containing a part of a bifurcation diagram
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

function Hampl = calcbranchHBM(sys,Section,Icalc,Hcalc)

Nb = length(Icalc);
nb_pt = size(Section.Upp,2);

st1 = Section.drawtype_init;
st2 = Section.drawtype_end;

H = sys.H;
DHp1 = 2*H+1;

if ~strcmp(st1,st2)
    ind_change = Section.ind_change-1;
    Hampl = zeros(Nb,nb_pt+1);
else
    Hampl = zeros(Nb,nb_pt);
end

for ii=1:Nb
    icalc = Icalc(ii);
    hcalc = Hcalc(ii);
    
    if icalc <= sys.nz
        I0   = (icalc-1)*DHp1+1;
        Icos = (icalc-1)*DHp1+1+hcalc;
        Isin = (icalc-1)*DHp1+1+H+hcalc;
    else
        I0   = sys.neq+1+ (icalc-1-sys.nz)*DHp1+1;
        Icos = sys.neq+1+ (icalc-1-sys.nz)*DHp1+1+hcalc;
        Isin = sys.neq+1+ (icalc-1-sys.nz)*DHp1+1+H+hcalc;
    end
    
    if hcalc == 0
        Hamplii = Section.Upp(I0,:);
    else
        Hamplii = sqrt(Section.Upp(Icos,:).^2+Section.Upp(Isin,:).^2);
    end
    
    if strcmp(st1,st2)
        Hampl(ii,:) = Hamplii;
    else
        Hampl_stab   = sqrt(Section.Ustab(Icos,:).^2+Section.Ustab(Isin,:).^2);
        Hampl(ii,:) = [Hamplii(1:ind_change-1),Hampl_stab,Hamplii(ind_change:end)];
    end
end

end
