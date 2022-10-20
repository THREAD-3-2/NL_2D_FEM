function [X,Y,Z] = extractXYdata(filename)
%%To unplot data from a matlab figure (.fig) files generated
% using version 7 or later. It can be used for both 2D and 3D plots
% Example usage:
% To unplot 2D graphs 
% [x,y] = unplot('example2D.fig')
% To unplot 3D graphs 
% [x,y,z] = unplot('example3D.fig') 
% Pradyumna
% January 2012
if nargin==1
    fig1 = load (filename,'-mat');
    a = fig1.hgS_070000.children(1).children(1).properties;
    index = length(fig1.hgS_070000.children(1).children(1).properties.XData);
    j = 1; % odd indexing parameter
    for i = 1:floor(length(fig1.hgS_070000.children(1).children)/2)
        if fig1.hgS_070000.children(1).children(j).type(1) == 't'
            break
        end
        if isfield(a,'ZData')
            % unploting 3D plot
            Y(index*(i-1)+1:index*i) = fig1.hgS_070000.children(1).children(j).properties.YData;
            X(index*(i-1)+1:index*i) = fig1.hgS_070000.children(1).children(j).properties.XData;
            Z(index*(i-1)+1:index*i) = fig1.hgS_070000.children(1).children(j).properties.ZData;
        else
            % unploting 2D plot
            Y(index*(i-1)+1:index*i) = fig1.hgS_070000.children(1).children(j).properties.YData;
            X(index*(i-1)+1:index*i) = fig1.hgS_070000.children(1).children(j).properties.XData;
        end
        j = j + 2;
    end 
else
    disp('Usage unplot(''filename.fig''). See ''help unplot'' for more details');
end