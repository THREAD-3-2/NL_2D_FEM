function [U0, omega0, lambda0,  qs_full, freq, shape, idx, amp] = initialise_MAN_computation(obj, sys, type, target_mode, varargin)
% Automatically initialise a MAN computation depending on the provided type
% type = 'autonomous' or 'forced'
%
%
H = sys.H;
% static solution
[qs_full, res] = obj.solve_static_problem();
% modal analysis around static configuration
[shape, freq] = obj.linear_modal_analysis(qs_full);
if strcmp(type,'forced')
    % deal with variable arguments
    if ~isempty(varargin)
        omega0 =varargin{1};
    else
        omega0 = freq(target_mode)*0.8*2*pi;
    end
    % Compute inital point for FRF (fixed frequency)
    lambda0 = omega0; % continuation parameter initial value
    [qp_full, bode] = obj.linear_analysis(H, omega0);
    [Z0] = obj.man_initial_point(H, omega0, qs_full, -qp_full);
    U0 = sys.init_U0(Z0, omega0, lambda0);
    U0 = obj.solve_MAN_system_at_fixed_frequency(U0, omega0, sys);
    % for outputs only (same as in 'autonomous')
    [val, idx] = max(shape(obj.boundary.active_dof,target_mode));
    idx = (2*H+1)*(idx-1) + 2; % take the first harmonic of dof that has the highest displacement amplitude in the shape
    amp = 1e-4*max( max(obj.mesh.nodes(:,2:3))- min(obj.mesh.nodes(:,2:3)) ); % 1e-4 * charactaritic length
elseif strcmp(type, 'autonomous')
    % deal with variable arguments
    if ~isempty(varargin)
        idx = varargin{1};
        amp = varargin{2};
    else % default parameter for NNM initialisation
        [val, idx] = max(shape(obj.boundary.active_dof,target_mode));
        idx = (2*H+1)*(idx-1) + 2; % take the first harmonic dof that has the highest amplitude in the shape
        amp = 1e-4*max( max(obj.mesh.nodes(:,2:3))- min(obj.mesh.nodes(:,2:3)) ); % 1e-4 * charactaritic length
    end
    % Compute inital point for NNM (fixed amplitude)
    omega0 = freq(target_mode)*2*pi;
    lambda0 = 0; % continuation parameter initial value
    epsilon = 1e-5; % amplitude of the linear mode
    [Z0] = obj.man_initial_point(H, omega0, qs_full, epsilon*shape(:,target_mode));
    U0 = sys.init_U0(Z0, omega0, lambda0);
    U0 = obj.solve_MAN_system_at_fixed_amplitude(U0, idx, amp, sys);
end



end
