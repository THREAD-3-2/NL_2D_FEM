function [ obj, params ] = displayvariables( obj, params )


nbmaxd = 8;
prompt = {'x1,y1','x2,y2','x3,y3','x4,y4','x5,y5','x6,y6','x7,y7','x8,y8'};
def =    {'0,0' ,'0,0' ,'0,0' ,'0,0' ,'0,0' ,'0,0' ,'0,0' ,'0,0' ,'0,0'};

%  user change dispars using the box
dispvars=params.dispvars;

for i=1:length(dispvars(:,1));
  def{i} = [num2str(dispvars(i,1)) ',' num2str(dispvars(i,2))];
end

rep = inputdlg(prompt,'Display variables',1,def);

dispvars = [];
for i=1:nbmaxd
  a = sscanf(rep{i},'%i,%i');
  if (a(1) <=0)
    break;
  end
  dispvars(i,1) = a(1); dispvars(i,2) = a(2);
end

params.dispvars=dispvars;

% recompute data for the proj bif diag and redraw each checkPoint

figure(2); clf; hold on;
for i=1:length(obj.ChckPoint)
    obj.ChckPoint{i}.dispvars = dispvars;
    obj.ChckPoint{i}.ncurve   = length(dispvars(:,1));
    obj.ChckPoint{i}.X = obj.ChckPoint{i}.Upp(dispvars(:,1),:)';   % X coordinate of the point to be plotted on Proj Bif Diag
    obj.ChckPoint{i}.Y = obj.ChckPoint{i}.Upp(dispvars(:,2),:)';   % X coordinate of the point to be plotted on Proj Bif Diag
    
    switch obj.ChckPoint{i}.BifStatus{1}
        case 'simplebif'
            obj.ChckPoint{i}.Xbif=obj.ChckPoint{i}.Ubif(dispvars(:,1));
            obj.ChckPoint{i}.Ybif=obj.ChckPoint{i}.Ubif(dispvars(:,2));
        otherwise
            if ~strcmp('nothing',obj.ChckPoint{i}.BifStatus{1})
                obj.ChckPoint{i}.Xstab=obj.ChckPoint{i}.Ustab(dispvars(:,1));
                obj.ChckPoint{i}.Ystab=obj.ChckPoint{i}.Ustab(dispvars(:,2));
            end
            
    end
    disp(obj.ChckPoint{i});
end