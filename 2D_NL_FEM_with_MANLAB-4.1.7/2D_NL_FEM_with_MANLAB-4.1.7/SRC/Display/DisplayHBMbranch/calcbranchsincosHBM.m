% calcbranchHBM   Branch computation for Manlab/HBM.
% function [Ampcos,Ampsin] = calcbranchsincosHBM(sys,Section,Icalc,Hcalc)
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
%    It returns [Ampcos,Ampsin], matrix of size (N,nb_pts),
%    containing the values of the branch sin and cos amplitudes. 
%    where nb_pts is the number of points to plot the projected
%    bifurcation diagram.
%
%    If the amplitude of the zeroth harmonics is asked, it return it in
%    Ampcos, with corresponding value of Ampsin equal to NaN.
%
%    By L. Guillot and O. Thomas / Nov. 2018

function [Ampcos,Ampsin] = calcbranchsincosHBM(sys,Section,Icalc,Hcalc)

Nb = length(Icalc);
nb_pt = size(Section.Upp,2);

st1 = Section.drawtype_init;
st2 = Section.drawtype_end;

H = sys.H;
DHp1 = 2*H+1;

if ~strcmp(st1,st2)
    ind_change = Section.ind_change-1;
    Ampcos = zeros(Nb,nb_pt+1);
    Ampsin = zeros(Nb,nb_pt+1);
else
    Ampcos = zeros(Nb,nb_pt);
    Ampsin = zeros(Nb,nb_pt);
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
        Ampcosii = Section.Upp(I0,:);
        Ampsinii = NaN*Section.Upp(I0,:);
    else
        Ampcosii = Section.Upp(Icos,:);
        Ampsinii = Section.Upp(Isin,:);
    end
    
    if strcmp(st1,st2)
        Ampcos(ii,:) = Ampcosii;
        Ampsin(ii,:) = Ampsinii;
    else
        Ampcos_stab   = Section.Ustab(Icos,:);
        Ampsin_stab   = Section.Ustab(Isin,:);
        Ampcos(ii,:) = [Ampcosii(1:ind_change-1),Ampcos_stab,Ampcosii(ind_change:end)];
        Ampsin(ii,:) = [Ampsinii(1:ind_change-1),Ampsin_stab,Ampsinii(ind_change:end)];
    end
end

end
