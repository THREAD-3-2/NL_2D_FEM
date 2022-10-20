function [val] = get(obj,propName)
% Get obj properties from the specified object and return the value
switch propName
    case 'U0now'
        val = obj.CurPoint.U0now;
    case 'Utnow'
        val = obj.CurPoint.Utnow; 
    case 'Ut2'
        val = obj.CurPoint.Ut2;   
    case 'Status'
        val = obj.CurPoint.Status;
    case 'ChckPoint'
        val = obj.ChckPoint;
    
end
