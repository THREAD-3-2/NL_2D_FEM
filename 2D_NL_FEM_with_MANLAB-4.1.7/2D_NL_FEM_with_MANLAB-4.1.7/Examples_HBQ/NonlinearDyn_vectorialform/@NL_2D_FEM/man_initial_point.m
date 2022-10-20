function [Z0] = man_initial_point(obj, H, omega0, qs_full, qp_full)
% compute MAN vector of variables from a given static configuration
% (qs_full) ad a given periodic config (qp_full).

qs =qs_full(obj.boundary.active_dof);
qp =qp_full(obj.boundary.active_dof,:);

[auxs] = obj.man_auxiliary_variables_vector(qs_full);

% [auxp] = obj.man_auxiliary_variables_vector(qp_full);
% [strain0, stress0] = obj.strains_and_stress_at_gauss_point(qs_full);
% [strainp, stressp] = obj.strains_and_stress_at_gauss_point(qp_full);

nz = 2*length(obj.boundary.active_dof);
nzaux = 12*obj.mesh.number_elements;
nz_tot = nz+nzaux;
number_elements = obj.mesh.number_elements;

Z0 = zeros(2*H+1, nz_tot);
% static part
Z0(1, 1:nz/2) = qs; % initial position
Z0(1, nz + (1:12*number_elements)) = auxs; % equations_vector.m
% cosine
for hh=1:size(qp,2)
    Z0(hh+1,1:nz/2) = real(qp(:,hh))'; % initial position
    Z0(hh+1+H,1:nz/2) = -imag(qp(:,hh))'; % initial position
    Z0(hh+1+H,nz/2+1:nz) = -omega0*real(qp(:,hh))'; % initial position
    Z0(hh+1,nz/2+1:nz) = omega0*imag(qp(:,hh))'; % initial position
%    Z0(hh+1, nz + (1:4*number_elements)) = auxp(1:4*number_elements)  ; % up, wp, thp theta 
%    Z0(hh+1+H, nz + (1:4*number_elements)) = -omega0*auxp(1:4*number_elements)  ; % up, wp, thp theta 
end
% Z0(2,nz + (4*number_elements + 1:5*number_elements))   = 0;      % theta at gauss point
% Z0(2,nz + (5*number_elements + 1:6*number_elements))   = strainp.meanrot;      % theta at gauss point
% Z0(2,nz + (6*number_elements + 1:7*number_elements))   = strainp.grad.up;    % theta at gauss point
% Z0(2,nz + (7*number_elements + 1:8*number_elements))   = strainp.grad.wp-strainp.meanrot;    % theta at gauss point
% Z0(2,nz + (8*number_elements + 1:9*number_elements))   = E*A*strainp.grad.up;     % theta at gauss point
% Z0(2,nz + (9*number_elements + 1:10*number_elements))  = k*G*A*(strainp.grad.wp-strainp.meanrot);     % theta at gauss point
% Z0(2,nz + (10*number_elements + 1:11*number_elements)) = E*I*strainp.grad.thetap;      % theta at gauss point
% Z0(2,nz + (11*number_elements + 1:12*number_elements)) = k*G*A*(strainp.grad.wp-strainp.meanrot);     % theta at gauss point
% 
% Z0(2+H,nz + (4*number_elements + 1:5*number_elements))   = 0;      % theta at gauss point
% Z0(2+H,nz + (5*number_elements + 1:6*number_elements))   = -omega0*strainp.meanrot;      % theta at gauss point
% Z0(2+H,nz + (6*number_elements + 1:7*number_elements))   = -omega0*strainp.grad.up;    % theta at gauss point
% Z0(2+H,nz + (7*number_elements + 1:8*number_elements))   = -omega0*strainp.grad.wp-strainp.meanrot;    % theta at gauss point
% Z0(2+H,nz + (8*number_elements + 1:9*number_elements))   = -omega0*E*A*strainp.grad.up;     % theta at gauss point
% Z0(2+H,nz + (9*number_elements + 1:10*number_elements))  = -omega0*k*G*A*(strainp.grad.wp-strainp.meanrot);     % theta at gauss point
% Z0(2+H,nz + (10*number_elements + 1:11*number_elements)) = -omega0*E*I*strainp.grad.thetap;      % theta at gauss point
% Z0(2+H,nz + (11*number_elements + 1:12*number_elements)) = -omega0*k*G*A*(strainp.grad.wp-strainp.meanrot);     % theta at gauss point

end