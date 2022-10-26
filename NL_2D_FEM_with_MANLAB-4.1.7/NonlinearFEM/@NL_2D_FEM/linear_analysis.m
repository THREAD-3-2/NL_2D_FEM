function [qp_full, bode] = linear_analysis(obj, H, Omega, varargin)
% Linear forced response analysis computes the amplitude and phase of the dof for each harmonics in the
% forcing vectors

% retreive assembled force vectors
fs_full = obj.vectors.static_forces;
fp_full = obj.vectors.periodic_forces;
% retreive assembled matrices
M_full = obj.matrices.mass;
C_full = obj.matrices.damping;
if ~isempty(varargin)
    qs_full = varargin{1}; % static solution
    K_full = obj.assemble_constant_matrix('stiffness_at_qs', qs_full); % tangent stiffness matrix
else
    K_full = obj.assemble_constant_matrix('stiffness_at_origin'); % stiffness matrix at origin
end
% apply bc
active_dof = obj.boundary.active_dof;
M = M_full(active_dof, active_dof);
C = C_full(active_dof, active_dof);
K = K_full(active_dof, active_dof);
fs = fs_full(active_dof);
fp = fp_full(active_dof,:);

Nw = length(Omega);
% compute linear solution (complex amplitude)
if Nw ==1 % if computation at a single frequency omega, return a matrix (size ndof x H)
    qp_full = zeros(3*obj.mesh.number_nodes, H);
    for h=1:size(fp,2)        
        w = Omega;
        Z = K-(h*w).^2*M+1i*w*h*C;
        X = Z\fp(:,h);
        qp_full(active_dof,h) = X;
        bode.amp_qp_full(active_dof, h) = abs(X);
        bode.phase_qp_full(active_dof, h) = angle(X);
        bode.qcos(active_dof, h) = real(X);
        bode.qsin(active_dof, h) = - imag(X);
    end
else % if computation over a frequency range, returns a cell array (length H) of matrices (size ndof x nW)
    qp_full = cell(H,1);
    for h=1:size(fp,2)
        qp_full{h} = zeros(3*obj.mesh.number_nodes, Nw);
        for ii=1:Nw
            w = Omega(ii);
            Z = K-(h*w).^2*M+1i*w*h*C;
            X = Z\fp(:,h);
            qp_full{h}(active_dof,ii) = X;
            bode.amp_qp_full{h}(active_dof,ii) = abs(Z\fp(:,h));
            bode.phase_qp_full{h}(active_dof,ii) = angle(Z\fp(:,h));
            bode.qcos{h}(active_dof, ii) = real(X);
            bode.qsin{h}(active_dof, ii) = -imag(X);
        end
    end
end

if 0 % TODO ploting functions
figure
subplot(2,1,1)
hold on
plot(Omega, bode.amp_qp_full{1}(4,:)) % u
plot(Omega, bode.amp_qp_full{1}(5,:)) % v
plot(Omega, bode.amp_qp_full{1}(6,:)) % theta
xlabel('Omega'); ylabel('Amp H1')
subplot(2,1,2)
hold on
plot(Omega, bode.phase_qp_full{1}(4,:)) % u
plot(Omega, bode.phase_qp_full{1}(5,:)) % v
plot(Omega, bode.phase_qp_full{1}(6,:)) % theta
xlabel('Omega'); ylabel('Phase H1')
end
end