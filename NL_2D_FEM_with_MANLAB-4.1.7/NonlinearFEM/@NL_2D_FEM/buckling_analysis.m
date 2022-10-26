function [frq_list, Lambda_range, N_list, q_list] = buckling_analysis(obj, varargin)
% Very simple buckling analysis. One computes the static equilibrium  for a
% static force of the type lambda*F_static (forces multiplies by factor lambda).
% Then a modal analysis is carried out and the smallest eigen value is computed.
% The list of lambda can be given or varargin. The default lambda list is of the type
% lambda = [1 .. 2, .. 5]
% The process stops if a zero eigenvalue appears

% varargin: Lambda_list
if ~isempty(varargin)
    Lambda_range = varargin{1};
else % increase from static vector fint = lambda*F_static
    Lambda_range = (1+[linspace(0.1,1,5), linspace(1,4,10)]);
end

q0_full = obj.vectors.null_vector;
fprintf(1,'Buckling Analysis \n')
lambda_prev = 1;
for i=1:length(Lambda_range)
    lambda = Lambda_range(i);
    
    fprintf(1,'  - step : %u / %u , ',i,length(Lambda_range))
    fprintf(1,'lambda : %2.2e \n, ',lambda)
    
    [qs_full, res] = obj.solve_static_problem(q0_full, lambda, lambda_prev);
    [shape, freq, imagfreq] = obj.linear_modal_analysis(qs_full);    
    [strain, stress] = obj.strains_and_stress_at_gauss_point(qs_full);
    
    N_list(:,i) = stress.N;
    frq_list(:,i) = freq;
    q_list(:,i)= qs_full;
    q0_full = qs_full;
    lambda_prev = lambda;
    
    if norm(imag(imagfreq))>1e-4
        fprintf(1,'Complex eigen value in buckling analysis, lambda= %2.2e, freq(1)=%2.2e',lambda, freq(1))
        warning('Complex eigen value in buckling analysis')
        break
    end
end
fig=figure(95); fig.Name='First Frequency vs Load coeff';
hold on
[Lambda_range, idx] = sort(Lambda_range(1:i),'ascend');
frq_list = frq_list(:,1:i);
frq_list = frq_list(:,idx);
plot(Lambda_range, frq_list(1,:));
% plot(Lambda_range, frq_list(2,:));
xlabel('Load coefficient \lambda')
ylabel('First eigen-frequency [S^{-1}]')
grid on
drawnow
end