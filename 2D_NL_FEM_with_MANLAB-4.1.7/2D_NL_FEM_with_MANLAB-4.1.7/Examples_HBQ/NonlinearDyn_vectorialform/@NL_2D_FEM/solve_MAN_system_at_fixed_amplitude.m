function U02 = solve_MAN_system_at_fixed_amplitude(obj, U0, idx, amp, sys)
% solve system at stating amplitude
fprintf(1, 'Solving MAN system at fixed amplitude, idx:%u, amp:%2.2e \n ',idx,amp)
% solver used (fsolve or home made newton)
solver = obj.solver; % 'homemade', 'fsolve'

% fsolve option
% options = optimoptions('fsolve','SpecifyObjectiveGradient',true,...
%       'display','off',...
%       'maxiter', solver.n_iter_max,...
%       'FunctionTolerance', solver.tol_res,...
%       'StepTolerance', solver.tol_step,...
%       'Algorithm','levenberg-marquardt');% trust-region-dogleg

ok_res=1;
niter=0;
Nmax=10; % number of restart
while ok_res
    niter=niter+1;
    if strcmp(solver.type,'homemade')
        [U02, res2, flag2, Jac] = obj.newton_solver('MAN_HBM_fixed_amp', U0, solver, idx, amp, sys);
    elseif strcmp(solver.type,'fsolve')        
        [U02, res2, flag2] = fsolve(@(U)man_residual_fixed_amplitude(U, idx, amp, sys), U0, options);
    end
    fprintf(1, '    res: %2.2e, flag: %u \n', norm(res2), flag2)
    if flag2>0
        ok_res=0;
        fprintf(1,'HBM system at initial amplitude solved, residual %2.2e \n', norm(res2))
    else
        U0=U02;
        fprintf(1,'    restarting, ')
    end
    if niter>Nmax
        warning('couldnt solve initial systemn too many restarts (10)')
        break
    end
    
end