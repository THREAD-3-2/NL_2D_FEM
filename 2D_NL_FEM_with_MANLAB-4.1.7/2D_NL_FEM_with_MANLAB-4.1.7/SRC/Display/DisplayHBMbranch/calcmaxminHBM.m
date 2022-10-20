% calcmaxminHBM   Computes max and min values of periodic
%                 oscillations for Manlab/HBM.
%    [Max Min]=plotbranchHBM(sys,Section,Icalc) computes the
%    maximum and miniumum amplitude of the periodic oscillations
%    computed by the Manlab/HBM method.
%
%    sys is an object containing the informations about the system solved.
%
%    Section is an object containing a part of a bifurcation diagram
%    computed with Manlab.
%
%    Icalc contains the indice of the entry of the initial physical
%    unknown vector z(t) whos max and min values are to be computed
%
%    It returns in Max and Min row vectors the max and min values
%    of u(t) as a function of the angular frequency.
%    length(Max)=length(Min)=size(U,2)
%
%    Example: [Max,Min]=plotbranchHBM(H,Neq,U,1)
%             computes max and min values of z1.
%
%    By L. Guillot and O. Thomas / Nov. 2018

function [Max,Min]=calcmaxminHBM(sys,Section,Icalc)


Nb = length(Icalc);
nb_pt=size(Section.Upp,2);

Nech=100*sys.H;
Nech = 1e3;
timegrid=linspace(0,1,Nech);

st1 = Section.drawtype_init;
st2 = Section.drawtype_end;

if strcmp(st1,st2)
    Max=zeros(Nb,nb_pt);
    Min=zeros(Nb,nb_pt);
    
    for k=1:nb_pt
        Utime=calcperiodHBM(sys,Section.Upp(:,k),Icalc,timegrid);
        Max(:,k)=max(Utime);
        Min(:,k)=min(Utime);
    end
    
else
    ind_change = Section.ind_change;
    Max=zeros(Nb,nb_pt+1);
    Min=zeros(Nb,nb_pt+1);
    
    kind = 1;
    for k=1:nb_pt+1
        if k == ind_change
            Uk = Section.Ustab;
        else
            Uk = Section.Upp(:,kind);
            kind = kind+1;
        end
        Utime=calcperiodHBM(sys,Uk,Icalc,timegrid);
        Max(:,k)=max(Utime);
        Min(:,k)=min(Utime);
    end
    
end



