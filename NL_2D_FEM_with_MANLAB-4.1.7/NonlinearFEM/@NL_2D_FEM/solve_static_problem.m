function [qs_full, res] = solve_static_problem(obj,  varargin)
% Solve the satic problem
% varargin = q0_full, lambda, lambda_prev
nvarargin = length(varargin);
if ~isempty(varargin) % initial guess provided by user
    q0_full = varargin{1};
else % default initial guess
    q0_full = zeros(3*obj.mesh.number_nodes,1);
end
if nvarargin>1 % variable coeff for buckling
    lambda = varargin{2}; % current force coeff : fint(q) = labmda*F_static
    if ~isempty(varargin{3})
        lambda_previous = varargin{3}; % previous force coeff (is sequential continuation)
    else
        lambda_previous = 1;
    end
    kr = lambda_previous/lambda;
    ok_lambda = 1;
else % default force coeff to 1
    lambda = 1;
    ok_lambda = 0; % no buckling analysis
end
% solver used (fsolve or home made newton)
solver = obj.solver; % 'homemade', 'fsolve'

% fsolve option
% options = optimoptions('fsolve','SpecifyObjectiveGradient',true,...
%       'display','off',...
%       'maxiter', solver.n_iter_max,...
%       'FunctionTolerance', solver.tol_res,...
%       'StepTolerance', solver.tol_step,...
%       'Algorithm','levenberg-marquardt');% trust-region-dogleg

fprintf(1,'Computing static solution using fsolve \n')
q0 = q0_full(obj.boundary.active_dof); % discard BC
Nramp_max = 3; % TODO
n = 0;
fprintf(1,'  ramp over static force amplitude \n')
for ramp = linspace(0.1, 1, Nramp_max) % ramp over fore amplitude
    n = n+1;
    res0 = norm(static_residual(obj,q0,ramp));
    fprintf(1,'    - iteration : %u / %u , ',n,Nramp_max)
    fprintf(1,'ramp : %2.2e, \n ',norm(ramp))
    fprintf(1,'initial res : %2.2e, ',norm(res0))
    OK_res = 1;
    nr = 0; % number of restarted fsolve
    while OK_res
        if ok_lambda==0
            if strcmp(solver.type,'homemade')
                [qs, res, flag, Jac] = obj.newton_solver('static', q0, solver, ramp);
            elseif strcmp(solver.type,'fsolve')
                [qs, res, flag, output, Jac] = fsolve(@(q)static_residual(obj,q, ramp),q0, options);
            end
        else
            if strcmp(solver.type,'homemade')
                [qs, res, flag, Jac] = obj.newton_solver('static', q0, solver, lambda*(kr*(1-ramp)+ramp));
            elseif strcmp(solver.type,'fsolve')
                [qs, res, flag, output, Jac] = fsolve(@(q)static_residual(obj,q, lambda*(kr*(1-ramp)+ramp)),q0, options);
            end
        end
        fprintf(1,'final res : %2.2e, ',norm(res))
        fprintf(1,'fsolve flag : %u \n', flag)
        
        if flag>0
            OK_res=0;
        else % restart
            nr = nr+1;
            fprintf(1,'  [static resolution (iter %u / %u)] : restarting (%u) :',n, Nramp_max, nr)
            fprintf(1,'initial res : %2.2e, ',norm(res))
            q0=qs;
        end
        if nr>500
            fprintf(1,'\n')
            error('Could not solve static problem, try reducing force amplitude or increase number of ramp steps')
        end
    end
    q0 = qs;
end
qs_full = zeros(3*obj.mesh.number_nodes,1);
qs_full(obj.boundary.active_dof)= qs;
% TODO update model vectors:
% if ok_lambda==0
%     obj.set_static_solution(qs_full)
% end
end