function [] = disp( obj )
figure(2); clf; hold on;
for i=1:length(obj.ChckPoint);
    disp(obj.ChckPoint{i});
end
