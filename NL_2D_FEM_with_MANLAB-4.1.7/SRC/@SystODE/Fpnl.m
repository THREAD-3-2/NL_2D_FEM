function [Fpnl_tot] = Fpnl(sys,U,p)
% Compute  the r.h.s Fpnl at order p from U (Taylor class)

switch sys.type
    case 'Equ'
        Fpnl_tot=sys.Fpnl_Equ(U,p);
    case 'HBM'
        Fpnl_tot=sys.Fpnl_HBM(U,p);
    case 'QPHBM'
        Fpnl_tot=sys.Fpnl_QPHBM(U,p);
end
