function obj = CheckPointCancel(obj,prop,nb)
% cancel obj checkpoint (1 or all) and return the updated object
figure(2);
switch prop
    case 'all'
        if isempty(obj.ChckPoint{1}) == 0
            obj.ChckPoint = cell(0);
        end
    case 'last'
        if isempty(obj.ChckPoint{end-nb+1}) == 0
            for i=1:nb
            obj.ChckPoint(end) = [];
            end
        end
    case 'one'
        [whatfound,idxsection,Apoint] = getindexes(obj);
        
      switch whatfound
      case 'regular'  
        disp(['cancel CheckPoint number ' num2str(idxsection)]);
        s=length(obj.ChckPoint);
        
        if s > 1
            for i=1:idxsection-1
                if i==1
                    liste=cell(0); liste{1}=obj.ChckPoint{1};
                else
                    liste{i}=obj.ChckPoint{i};
                end
            end
            for i=idxsection:s(1)-1
                liste{i}=obj.ChckPoint{i+1};
            end
            obj.ChckPoint=liste;
        else
            obj.ChckPoint = cell(0);
        end
      end
end
