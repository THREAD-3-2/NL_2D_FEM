function [obj] = set(obj,varargin)
%  Set obj properties and return the updated object
propertyArgin = varargin;
while length(propertyArgin) >=2
    prop = propertyArgin{1};
    val=propertyArgin{2};
    propertyArgin=propertyArgin(3:end);
    switch prop
        case 'U0'
            obj.U0 = val;
        case 'Ut'
            obj.Ut = val;
        case 'Amax'
            obj.Amax = val;
        case 'Uend'
            obj.Uend = val;
        case 'Utend'
            obj.Utend = val;
        case 'Upp'
            obj.Upp = val;
        case 'BifStatus'
            obj.BifStatus = val;
        case 'Ubif'
            obj.Ubif = val;
        case 'Utf'
            obj.Utf = val;
        case 'Utb'
            obj.Utb = val;
        case 'alpha'
            obj.alpha = val;
        case 'Uscale'
            obj.Uscale = val;
        case 'dispcolors'
            obj.dispcolors = val;
        case 'dispvars'
            obj.dispvars = val;
        case 'drawtype'
            obj.dratype = val;
        case 'markers'
            obj.markers = val;
        case 'nbpts'
            obj.nbpts = val;
        case 'ncurve'
            obj.ncurve = val;
        case 'X'
            obj.X = val;
        case 'Y'
            obj.Y = val;
        case 'A'
            obj.A = val;
        case 'Xbif'
            obj.Xbif = val;
        case 'Ybif'
            obj.Ybif = val;
        case 'Xstab'
            obj.Xstab = val;
        case 'Ystab'
            obj.Ystab = val;
        case 'Achange'
            obj.Achange = val;
        case 'Eigen'
            obj.Eigen = val;
    end
end


