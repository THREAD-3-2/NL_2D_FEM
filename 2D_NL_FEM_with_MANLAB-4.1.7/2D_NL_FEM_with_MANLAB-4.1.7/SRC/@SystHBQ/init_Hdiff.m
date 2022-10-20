function [Unew] = init_Hdiff(sys,Uold)

Hold = ((size(Uold,1)-4)/sys.nz_tot-1)/2;

neq_old = sys.nz*(2*Hold+1);
neq_aux_old = sys.nz_aux*(2*Hold+1);
Zold = reshape(Uold([1:neq_old neq_old+3:neq_old+3+neq_aux_old-1]) , 2*Hold+1 , sys.nz_tot); 
omega = Uold(neq_old+1);
lambda = Uold(neq_old+2);

Hnew = sys.H;
DHp1 = 2*Hnew+1;

if Hold >= Hnew
    Znew = Zold([(1:Hnew+1) (Hold+2:Hold+1+Hnew)],:);
else
    Znew = zeros(DHp1,sys.nz_tot);
    Znew([(1:Hold+1) (Hnew+2:Hnew+1+Hold)],:) = Zold;
end

Unew = sys.init_U0(Znew,omega,lambda);

end

