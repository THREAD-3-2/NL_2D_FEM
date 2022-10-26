function [U0] = init_U0(sys,Z,omega1,omega2,lambda)

H = sys.H;
nb_coef = 2*((2*H+1)*H + H)+1;

Ufourier = reshape(Z(:,1:sys.nz),nb_coef*sys.nz,1);
Uaux_fourier = reshape(Z(:,sys.nz+1:end),nb_coef*sys.nz_aux,1);

U0 = [Ufourier;omega1;omega2;lambda;Uaux_fourier;omega1^2;omega2^2;omega1*omega2;lambda*omega1;lambda*omega2];

end

