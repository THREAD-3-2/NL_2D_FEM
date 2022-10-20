% plotbarHBM   Spectrum of periodic solution for Manlab/HBM
%    Hand = plotbarsincosHBM(sys,U,Idisp,fig) plots the spectrum,
%    as a stem like graph of the harmonics content of the solution
%    computed by the Manlab/HBM method.
%
%    sys is an object containing the informations about the system solved.
%
%    U is the vector of unknowns of the final algebraic system solved
%    by Manlab.
%
%    Idisp contains the indices of the entries of the initial physical
%    unknown vector u(t) to be displayed.
%
%    It returns a vector of barseries handles in Hand. Each entry of
%    Hand correspond to a bar
%
%    Example: Hand=plotbarsincosHBM(sys,U,[1 2],1)
%
%    By L. Guillot and O. Thomas / Nov. 2018

function Hand=plotbarsincosHBM(sys,U,Idisp)

% Set of colors for display.
Col = dispcolorML;

% Column vector test
% ==================

if size(U,2)>size(U,1)
    U=U';
end

if size(U,2)>1
    U=U(:,1);
    disp('plotbarsincosU: U should be a column vector')
end

% Plotting the bargraph
% =====================

Hand=[];
Ndisp=length(Idisp);
H = sys.H;
DHp1 = 2*H+1;

for ii=1:Ndisp
    idisp=Idisp(ii);
    cdisp=Col(ii,:);
    
    if idisp <= sys.nz
        I0   = (idisp-1)*DHp1+1;
        Icos = (idisp-1)*DHp1+1+(1:H);
        Isin = (idisp-1)*DHp1+1+H+(1:H);
    else
        I0   = sys.neq+1+ (idisp-1-sys.nz)*DHp1+1;
        Icos = sys.neq+1+ (idisp-1-sys.nz)*DHp1+1+(1:H);
        Isin = sys.neq+1+ (idisp-1-sys.nz)*DHp1+1+H+(1:H);
    end
    
    Xtick=ii:Ndisp:(Ndisp*2*H+ii);
    Ubar= [U(Icos),U(Isin)]';
    Ubar=[U(I0);Ubar(:)];
    hand=bar(Xtick,Ubar);
    set(hand,'Facecolor',cdisp,'barwidth',.8/Ndisp)
    Hand=[Hand;hand];
    Leg{ii}=['z_',num2str(idisp)];
end

% Labels
% ======

Xtick=1:Ndisp:(2*H+1)*Ndisp;
Ntick=length(Xtick);
Xtickstr=cell(1,Ntick);
Xtickstr{1}='H0';
for ii=1:H
    Xtickstr{2*ii}=['H',num2str(ii),'c'];
    Xtickstr{2*ii+1}=['H',num2str(ii),'s'];
end

dI=2*ceil(H/6)-2;
if dI>1
    Xtickstr0=cell(1,Ntick);
    Xtickstr0(2:dI:end)=Xtickstr(2:dI:end);
    Xtickstr=Xtickstr0;
end

set(gca,'Xtick',Xtick,'XtickLabel',Xtickstr,'Xlim',[0 (2*H+1)*Ndisp+1],'Xgrid','on')
ylabel('sin/cos amplitude')
box on

legend(Hand,Leg);



