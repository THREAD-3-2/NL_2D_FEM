function [Znew,omega1,omega2,lambda] = init_Hdiff(sys,Uold)

Hold = floor(sqrt(((size(Uold,1)-8)/sys.nz_tot-1)/4));
nb_coef_old = 2*((2*Hold+1)*Hold + Hold)+1;

neq_old = sys.nz*nb_coef_old;
neq_aux_old = sys.nz_aux*nb_coef_old;
Zold = reshape(Uold([1:neq_old neq_old+4:neq_old+4+neq_aux_old-1]) , nb_coef_old , sys.nz_tot); 
omega1 = Uold(neq_old+1);
omega2 = Uold(neq_old+2);
lambda = Uold(neq_old+3);

Hnew = sys.H;
nb_coef_new = 2*((2*Hnew+1)*Hnew + Hnew)+1;
    
if Hnew >= Hold
    Znew = zeros(nb_coef_new,sys.nz_tot);
    indices_cos = [(1:Hold+1)' ; reshape(repmat(2*Hnew+1 + 1 + (-Hold:1:Hold)',1,Hold) + (0:1:Hold-1)*(2*Hnew+1),(2*Hold+1)*Hold,1)];
    indices_sin = indices_cos(2:end)+(2*Hnew+1)*Hnew + Hnew;
    Znew([indices_cos;indices_sin],:) = Zold;
else
    disp('Pas encore implémenté.');
end

end

