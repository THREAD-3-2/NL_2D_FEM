function Section = CheckPointExport(obj,prop)
% replace the "current point" in section idxsection for a=Apoint
%For complex bifurcation point seting, an additional condition must be
%added.
global Ck
switch prop
    case 'all'
        Section = obj.ChckPoint;
    case 'one'
        [whatfound,idxsection,Apoint] = getindexes(obj);
        Section = obj.ChckPoint{idxsection};
end
end
