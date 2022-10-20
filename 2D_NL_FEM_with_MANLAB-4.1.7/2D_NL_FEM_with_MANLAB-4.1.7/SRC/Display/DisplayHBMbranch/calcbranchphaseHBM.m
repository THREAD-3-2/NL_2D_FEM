% calcbranchphaseHBM   Branch phase components computation for Manlab/HBM.
%    Phases=calcbranchphaseHBM(sys,Section,Icalc,Hcalc) computes the phase
%    harmonics components of the solution computed by the Manlab/HBM method.
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
%    It returns Phase, a matrix of size (length(Icalc),nb_pts),
%    where nb_pts is the number of points to plot the projected
%    bifurcation diagram.
%
%    By L. Guillot and O. Thomas / Nov. 2018

function Phase = calcbranchphaseHBM(sys,Section,Icalc,Hcalc)

Nb=length(Icalc);
nb_pt = size(Section.Upp,2);

st1 = Section.drawtype_init;
st2 = Section.drawtype_end;

H = sys.H;
DHp1 = 2*H+1;

if ~strcmp(st1,st2)
    ind_change = Section.ind_change-1;
    Phase = zeros(Nb,nb_pt+1);
else
    Phase = zeros(Nb,nb_pt);
end

for ii=1:Nb
    icalc = Icalc(ii);
    hcalc = Hcalc(ii);
    
    if icalc <= sys.nz
        Icos = (icalc-1)*DHp1+1+hcalc;
        Isin = (icalc-1)*DHp1+1+H+hcalc;
    else
        Icos = sys.neq+1+ (icalc-1-sys.nz)*DHp1+1+hcalc;
        Isin = sys.neq+1+ (icalc-1-sys.nz)*DHp1+1+H+hcalc;
    end
    
    if hcalc == 0
        Phaseii = NaN*ones(1,nb_pt);
    else
        Phaseii = -atan2(Section.Upp(Isin,:),Section.Upp(Icos,:));
    end
    
    if strcmp(st1,st2)
        Phase(ii,:) = Phaseii;
    else
        if hcalc == 0
            Phase(ii,:) = NaN*ones(1,nb_pt+1);
        else
            Phasestab = -atan2(Section.Ustab(Isin,:),Section.Ustab(Icos,:)) / pi;
            Phase(ii,:) = [Phaseii(1:ind_change-1),Phasestab,Phaseii(ind_change:end)];
        end
    end
end


end
