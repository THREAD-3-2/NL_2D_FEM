function [X, res, flag, Jac] = newton_solver(obj, type, X0, solver, varargin)
% Homemade newton solver
% X0 initial guess
% type = 'static', 'MAN_HBM_fixed_amp', 'MAN_HBM_fixed_freq'
% varargin according to type
% 'static': varargin empty or ramp coeff (1x1 double)
% 'MAN_HBM_fixed_amp', vararing = {idx, amp, sys}
% 'MAN_HBM_fixed_freq', varargin = {omega0, sys}
%
% flag, 1: solved at initial point, 2: solved after iteration, 0: number of iteration exceded    

%% set tolerances from solver parameters
tol_res = solver.tol_res;
tol_step = solver.tol_step;
n_iter_max = solver.n_iter_max;
fprintf(1,'system solving using homemade newton (tol_res : %2.2e, tol_step : %2.2e, n_iter_max : %u) \n', tol_res, tol_step, n_iter_max)

%% retreive residual function from type = 'static', 'MAN_HBM_fixed_amp', 'MAN_HBM_fixed_freq'
if strcmp(type, 'static')
    if isempty(varargin)
        func_res = @(q)obj.static_residual(q);
    else
        func_res = @(q)obj.static_residual(q, varargin{1}); % ramp coeff
    end
elseif strcmp(type, 'MAN_HBM_fixed_amp')
    if isempty(varargin)
        error('varargin should not be empty')
    else
        func_res = @(U0)man_residual_fixed_amplitude(U0, varargin{1}, varargin{2}, varargin{3} ); % idx, amp, sys
    end
elseif strcmp(type, 'MAN_HBM_fixed_freq')
    if isempty(varargin)
        error('varargin should not be empty')
    else
        func_res = @(U0)man_residual_fixed_frequency(U0, varargin{1}, varargin{2}); % omega0, sys
    end
end

%% initial residual and jacobian
[res, Jac] = func_res(X0);
norm_res = norm(res);
if norm_res < tol_res % if solved at initial point, return
    fprintf(1,'system solved at initial point, norm_res : %2.2e \n', norm_res)
    X=X0; flag = 1;
    return
end
fprintf(1,'    Initial residual, norm_res : %2.2e \n', norm_res)
norm_step = 1;

%% Newton iterations
n_iter = 0;
while (norm_res >= tol_res || norm_step >= tol_step) && n_iter <=n_iter_max
    n_iter = n_iter+1;
    % compute step
    delta_X = -Jac\res;
    X0 = X0 + 0.9*delta_X;
    % residual and jacobian at new point
    [res, Jac] = func_res(X0);
    norm_res = norm(res);
    norm_step = norm(delta_X);
    fprintf(1,'    iter %u, norm_step: %2.2e, norm_res : %2.2e \n', n_iter, norm_step, norm_res)
    
    % check for diverging residua TODO
    if n_iter > 1 && (norm_res/norm_resPrev) > 10
        warning(' increasing residual in home made fsolve \n')        
    end
    norm_resPrev = norm_res;
end

%% Outputs setings
if n_iter>=n_iter_max % too many iteration
    fprintf(1,'Error Too many newton iteration, couldnt solve problem (increase n_iter_max) \n')
    flag = 0;
    X = X0;
else % system solved
    fprintf(1, 'System solved \n')
    flag = 2;
    X = X0;
end
end
