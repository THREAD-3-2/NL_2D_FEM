function obj = CheckPointAdd(obj,ChckPoint)
% Add a new checkpoint
s = length(obj.ChckPoint);
obj.ChckPoint{s+1} = ChckPoint;

