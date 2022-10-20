function [U] = visupoint(obj,params)
global Ck

[whatfound,idxsection,Apoint] = getindexes(obj);
U0=get(obj.ChckPoint{idxsection},'U0');
U  = evalseries(U0,Apoint,Ck); 


