function plot_mode_shape(obj, mode_list, shapes, varargin)
% Plot a list of mode shapes
% varargin : fig, marker
marker = 'b-'; % default marker
% deals with variables argument
if ~isempty(varargin)
    fig=varargin{2}; % figure handle
    if length(varargin)>1
        marker=varargin{3}; % marker spec
    end
else
    fig = figure(98);
end
Nmode = length(mode_list);
shapes = shapes(:,mode_list);
fig.Name='Mode shapes'; hold on
Nlines = floor(Nmode/3)+1; % number of lines in the subplot matrix (3 columns)
for i=1:Nmode
    subplot(Nlines, 3, i);   hold on
    % plot undeformed mesh
    obj.plot_deformed_mesh(obj.vectors.null_vector, fig, '--k') % undef
    qs_full = shapes(:,i);
    obj.plot_deformed_mesh(qs_full, fig, marker) % def
    % set axes
    % find mesh limits
    maxXYmesh = max(obj.mesh.nodes(:,2:3));
    minXYmesh = min(obj.mesh.nodes(:,2:3));
    % find displacement limits
    maxXdisp = max(abs(qs_full(1:3:end)));
    maxYdisp = max(abs(qs_full(2:3:end)));
    % define axis limits
    axis([(minXYmesh(1)-maxXdisp)*1.1, (maxXYmesh(1)+maxXdisp)*1.1, (minXYmesh(2)-maxYdisp)*1.1, (maxXYmesh(2)+maxYdisp)*1.1])
end
end