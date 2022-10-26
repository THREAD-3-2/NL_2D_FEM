function plot_static_configuration(obj, qs_full, varargin)
% plot the static configuration
% varargin : fig, marker
marker = 'k-'; % default marker
% deals with variables argument
if ~isempty(varargin)
    fig=varargin{2}; % figure handle 
    if length(varargin)>1
        marker=varargin{3}; % marker spec
    end
else
    fig = figure;
end
% plot deformed mesh
fig.Name='Static Configuration'; hold on
obj.plot_deformed_mesh(obj.vectors.null_vector, fig, '--k') % undef
obj.plot_deformed_mesh(qs_full, fig, marker) % def
% set axes
% find mesh limits
maxXYmesh = max(obj.mesh.nodes(:,2:3));
minXYmesh = min(obj.mesh.nodes(:,2:3));
% find displacement limits
maxXdisp = max(abs(qs_full(1:3:end)));
maxYdisp = max(abs(qs_full(2:3:end)));
% define axis limits
if maxXdisp>0 && maxYdisp>0
axis([(minXYmesh(1)*1.1-maxXdisp)*1.1, (maxXYmesh(1)*1.1+maxXdisp)*1.1, (minXYmesh(2)*1.1-maxYdisp)*1.1, (maxXYmesh(2)*1.1+maxYdisp)*1.1])
end
end