function [obj] = set(obj,varargin)
% Set obj properties and return the updated object
propertyArgin = varargin;
while length(propertyArgin) >=2
    prop = propertyArgin{1};
    val=propertyArgin{2};
    propertyArgin=propertyArgin(3:end);
    switch prop
        case 'U0now'
            obj.CurPoint.U0now=val;
        case 'Utnow'
            obj.CurPoint.Utnow=val; 
        case 'Ut2'
            obj.CurPoint.Ut2=val;
        case 'Statut'
            obj.CurPoint.Status=val;
    end
end
