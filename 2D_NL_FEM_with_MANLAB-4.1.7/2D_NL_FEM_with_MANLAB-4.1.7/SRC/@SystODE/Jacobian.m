function Jacobian = Jacobian(sys,Utot)
% Compute the Jacobian of R with respect to U

switch sys.type
    case 'Equ'
        Jacobian=sys.Jacobian_Equ(Utot);
    case 'HBM'
        Jacobian=sys.Jacobian_HBM(Utot);
    case 'QPHBM'
        Jacobian=sys.Jacobian_QPHBM(Utot);
end
