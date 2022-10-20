function [shape, freq, imagfreq] = linear_modal_analysis(obj, varargin)
% linear modal analysis. can be arround reference config (varargin empty) 
% or around a given config (varargin{1} = qs).

active = obj.boundary.active_dof;
Nmode = length(active);
% FEM model linear matrices
M_full = obj.matrices.mass;
if isempty(varargin) % modal analysis with tangent stifness at origin
    K_full = assemble_constant_matrix(obj, 'stiffness_at_origin');
else % modal analysis with tangent stifness at a given state qs
    qs = varargin{1};
    K_full = assemble_constant_matrix(obj, 'stiffness_at_qs', qs);
end
% truncate to account for bc
K = K_full(active, active);
M = M_full(active, active);
% eigen value problem solving
[V, D] = eig(K,M);
test = diag(D)<0;
if sum(test)
    warning('there is negative eigenvalues, complex frequency')
%    return
end
imagfreq = imag(sqrt(diag(D)))/2/pi;
freq = real(sqrt(diag(D)))/2/pi;
% resizing mode shape to account for bc
shape = zeros(3*obj.mesh.number_nodes, Nmode);
shape(active,:)= V./repmat(max(abs(V)),size(V,1),1);
% shape(active,:)= V;
% sorting modes by ascending frequency
[freq, idx] = sort(freq, 'ascend');
shape = shape(:,idx);
end