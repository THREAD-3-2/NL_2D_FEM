% plotcoefQPHBM   Spectrum of quasi-periodic solution for Manlab/QPHBM 
%    Hand = plotcoefQPHBM(sys,U,Idisp) plots the spectrum,
%    of the harmonics content of the solution
%    computed by the Manlab/QPHBM method.
%
%    sys is an object containing the informations about the system solved.
%
%    U is the vector of unknowns of the final algebraic system solved
%    by Manlab. 
% 
%    Idisp contains the indices of the entries of the initial physical 
%    unknown vector u(t) to be displayed.
%
%    Example: figure; Hand=plotcoefQPHBM(sys,U,[1 2])
%
%    By L. Guillot and O. Thomas / Nov. 2018

function Hand=plotcoefQPHBM(sys,U,Idisp)

H = sys.H;

if numel(Idisp)>1
    disp('plotcoefQPHBM.m : Only the first variable is displayed. Sum of variables considered.');
end

Ztot = sys.get_Ztot(U);
Varcompl = sys.real_to_compl(sum(Ztot(:,Idisp),2));
Matcoef = sys.get_mat_var(Varcompl);

colormap jet;
Hand = imagesc(-H:H,-H:H,log10(abs(Matcoef)));          %           Peut-être flipud..
colorbar

end


