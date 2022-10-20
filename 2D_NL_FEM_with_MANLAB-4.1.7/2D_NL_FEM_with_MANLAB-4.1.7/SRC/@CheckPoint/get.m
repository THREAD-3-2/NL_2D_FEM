function [val] = get(obj,propName)
% GET Get branch properties from the specified object
% and return the value
switch propName
    case 'U0'
        val = obj.U0;
    case 'Ut'
        val = obj.Ut;
    case 'Amax'
        val = obj.Amax;
    case 'Uend'
        val = obj.Uend;
    case 'Utend'
        val = obj.Utend;
    case 'Upp'
        val = obj.Upp;
    case 'BifStatus'
        val=obj.BifStatus;
    case 'Ubif'
        val=obj.Ubif;
    case 'Utf'
        val=obj.Utf;
    case 'Utb'
        val=obj.Utb;
    case 'alpha'
        val=obj.alpha;
    case 'Uscale'
        val=obj.Uscale;
    case 'dispcolors'
        val = obj.dispcolors;
    case 'dispvars'
        val = obj.dispvars;
    case 'drawtype'
        val = obj.dratype;
    case 'markers'
        val = obj.markers;
    case 'nbpts'
        val = obj.nbpts;
    case 'ncurve'
        val = obj.ncurve;
    case 'X'
        val = obj.X;
    case 'Y'
        val = obj.Y;
    case 'A'
        val = obj.A;
    case 'Xbif'
        val=obj.Xbif;
    case 'Ybif'
        val=obj.Ybif;
    case 'Xstab'
        val = obj.Xstab;
    case 'Ystab'
        val = obj.Ystab;
    case 'Astab'
        val = obj.Astab;
    case 'Eigen'
        val = obj.Eigen;
end