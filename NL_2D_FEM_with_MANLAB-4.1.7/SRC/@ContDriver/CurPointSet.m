function [obj,params,sys] = CurPointSet(obj,params,sys)
% replace the "current point" in section idxsection for a=Apoint
%For complex bifurcation point seting, an additional condition must be
%added.
global Ck 

[whatfound,idxsection,Apoint] = getindexes(obj);

switch whatfound
    case 'simplebif'
        obj.CurPoint.U0now=get(obj.ChckPoint{idxsection},'Ubif');
        obj.CurPoint.Utnow=get(obj.ChckPoint{idxsection},'Utb');
        obj.CurPoint.Ut2 = get(obj.ChckPoint{idxsection},'Utf');
        obj.CurPoint.BifStatus='simplebif';
        
    otherwise
        U0 = get(obj.ChckPoint{idxsection},'U0');
        obj.CurPoint.U0now  = evalseries(U0,Apoint,Ck);
        obj.CurPoint.Utnow  = evalderiv(U0,Apoint,Ck);
        obj.CurPoint.BifStatus =whatfound;
end
end
