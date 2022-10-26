function [coord] = getcoord(sys,variablename,varargin)
% function [coord] = getcoord(sys,variablename,varargin)
% gives the coordinates of the variable 'variablename' in the total vector
% of unknowns.
%
% The third and fourth argument values specifies the number of the variable 
% and of the harmonics considered.
%
% Example : getcoord('lambda') would give coord = sys.neq+1.
%           getcoord('cos',1,5) would give the index of the coefficient =
% of the fifth cosine of the first variable : that is 6.

H = sys.H;
DHp1 = 2*H+1;

switch variablename
    case 'omega'
        coord = sys.neq;
    case 'lambda'
        coord = sys.neq+1;
    case 'cos'
        ivar = varargin{1};
        ih = varargin{2};
        if ivar <= sys.nz
            coord = (ivar-1)*DHp1 + ih+1;
        else
            coord = sys.neq + 1 + (ivar-sys.nz-1)*DHp1 + ih+1;
        end
    case 'sin'
        ivar = varargin{1};
        ih = varargin{2};
        if ivar <= sys.nz
            coord = (ivar-1)*DHp1 + H + ih+1;
        else
            coord = sys.neq + 1 + (ivar-sys.nz-1)*DHp1 + H + ih+1;
        end
end

