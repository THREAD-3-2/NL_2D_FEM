function x = horzcat(varargin)
%HORZCAT  of a Taylor object.
%
%   See also Taylor, Taylor/VERTCAT
%
%MATLAB Diamant toolbox.
%K. Lampoh,  Isabelle Charpentier B. Cochelin.

for counter=1:nargin
    if isa(varargin{counter}, 'Taylor')
        x.order=varargin{counter}.order;
    end
end

x.value=[];
x.coef=[];

for jj = 1:nargin
    if isa(varargin{jj}, 'Taylor')
        x.value=[x.value varargin{jj}.value];
        x.coef=[x.coef varargin{jj}.coef];
        x.order=varargin{jj}.order;

    elseif ~isempty(varargin{jj})
        x.value=[x.value varargin{jj}];
        x.coef=[x.coef zeros([size(varargin{jj}),x.order])];
    end
end

x=Taylor(x.order,x.value,x.coef);