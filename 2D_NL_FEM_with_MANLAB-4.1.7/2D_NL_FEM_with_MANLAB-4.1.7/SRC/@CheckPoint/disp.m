function disp(obj)
figure(2); hold on

for i=1:obj.ncurve        
    
    if obj.Astab > 0
        plot([obj.X(1:obj.ind_change-1,i) ; obj.Xstab(i)],[obj.Y(1:obj.ind_change-1,i) ; obj.Ystab(i)], strcat(obj.dispcolors(i), obj.drawtype_init), ...
            [obj.Xstab(i) ; obj.X(obj.ind_change:end,i)],[obj.Ystab(i) ; obj.Y(obj.ind_change:end,i)], strcat(obj.dispcolors(i), obj.drawtype_end)); 
        
%         plot([obj.X(1:obj.ind_change,i) ; obj.Xstab(i)],[obj.Y(1:obj.ind_change,i) ; obj.Ystab(i)], strcat(obj.dispcolors(i), obj.drawtype_init), ...
%             [obj.Xstab(i) ; obj.X(obj.ind_change+1:end,i)],[obj.Ystab(i) ; obj.Y(obj.ind_change+1:end,i)], strcat(obj.dispcolors(i), obj.drawtype_end)); 
    else
       plot(obj.X(:,i),obj.Y(:,i), strcat(obj.dispcolors(i), obj.drawtype_end)); 
    end
    
    if (obj.markers==1)
        plot(obj.X(1,i),obj.Y(1,i), strcat(obj.dispcolors(i),'.')); 
    end
    
    for i_bif=1:numel(obj.BifStatus)
        switch obj.BifStatus{i_bif}
            case 'simplebif'
                plot(obj.Xbif(i),obj.Ybif(i),strcat(obj.dispcolors(i),'o') );
            case 'nothing'
                disp('');
            otherwise
                plot(obj.Xstab(i),obj.Ystab(i),strcat(obj.dispcolors(i),'p') );
        end
    end
    
end
