function [U0] = init_U0(sys,Z,omega,lambda)

H = sys.H;
DHp1 = 2*H+1;

if numel(Z) == DHp1*sys.nz_tot
    Ufourier = reshape(Z(:,1:sys.nz),DHp1*sys.nz,1);
    Uaux_fourier = reshape(Z(:,sys.nz+1:end),DHp1*sys.nz_aux,1);
% else
%     Uauxt = zeros(DHp1,obj.nz_aux);
%     for tt=1:DHp1
%         t = 2*pi/omega*(tt-1)/DHp1;
%         VectCos = cos((1:H)*omega*t);
%         VectSin = sin((1:H)*omega*t);
%         MAT=[ 1 , VectCos , VectSin ];
%         Ut=MAT*Z;
%         dUt=MAT*omega*obj.D(Z);
%         Uauxt(tt,:) = obj.auxiliary_variables(obj,[Ut lambda]',dUt');
%     end
%     FFT_Uaux = fft(Uauxt)/(2*H+1);
%     Uaux = obj.compl_to_real(FFT_Uaux(1:H+1,:));
%     
%     Ufourier = reshape(Z,DHp1*obj.nz,1);
%     Uaux_fourier = reshape(Uaux,DHp1*obj.nz_aux,1);
%     
end


U0 = [Ufourier;omega;lambda;Uaux_fourier;omega^2;lambda*omega];

end

