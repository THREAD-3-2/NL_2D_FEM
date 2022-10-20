function [res, jac] = static_residual(obj, q, varargin)
% static residual: res = fint(q) - f_static

if ~isempty(varargin)
    ramp = varargin{1}; % coeff between 0 and 1
else
    ramp = 1;
end
% q is of length ndof - n_constrained_dof
active_dof = obj.boundary.active_dof;
q_full = zeros(3*obj.mesh.number_nodes, 1);
q_full(active_dof) = q;
% internal forces
fint = obj.assemble_constant_vector('internal_force_vector', q_full);
% residual
res = fint - ramp*obj.vectors.static_forces;
% apply boundary condition
res = res(obj.boundary.active_dof); 
% Jacobian computation if needed
if nargout>1
    jac = obj.assemble_constant_matrix('stiffness_at_qs', q_full);
    jac = jac(active_dof, active_dof);
end

end