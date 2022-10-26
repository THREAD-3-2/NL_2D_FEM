% plotNNMHBM   Nonlinear Normal Mode / Invariant Manifold plot
%    Hand = plotNNMHBM(sys,U,Idisp,Ia) plots the manifold to
%    which the periodic solutions computed by the Manlab/HBM method
%    belongs.
%
%    U is the vector of unknowns of the final algebraic system solved
%    by Manlab. 
% 
%    Idisp contains the indices of the entries of the initial physical 
%    unknown vector u(t) to be displayed. Since the manifold is
%    displayed in 3D, length(Idisp) must be 3. 
%    If Idisp=[1 3 2], z1 is the x-axis, z3 is the y-axis and z2 is
%    the zaxis. 
%
%    Ia contains the indices of the arclength of the branch of
%    periodic solutions to be displayed. It correspond to the
%    radial coordinate of the displayed manifold. 
%    Ia must be such that min(Ia)>=1 and max(Ia)<=size(U,2)
%    
%    fig is a figure/axis handle to plot the manifold on.
%
%    It returns in Hand the handle to the displayed surf object. 
%
%    Example: Hand=plotNNMHBM(Neq,H,U,[1 3 2],1:10:200,1)
%
%    By L. Guillot and O. Thomas / Nov. 2018


function Hand=plotNNMHBM(sys,Section,Idisp,Ia)

if length(Idisp)~=3
  disp('length(Idisp) should be 3')
  Hand=[];
  return
end 

time = linspace(0,1,100);
time = time(:);

Na=length(Ia);
Nech = length(time);
Ndisp = length(Idisp);

if isfloat(Section)
    U = Section;
else
    U = Section.Upp;
end

if min(Ia)>size(U,2)
  disp('min(Ia) should be <= size(U,2)')
  Ia=1:size(U,2);
end
 
if max(Ia)>size(U,2)
  disp('max(Ia) should be <= size(U,2)')
  Ia=1:size(U,2); 
end

% Calcul des orbites périodiques
% ==============================

Uplot=zeros(Nech,Na,Ndisp);

for ii=1:Ndisp
  idisp=Idisp(ii);
  for jj=1:Na
    ia=Ia(jj);
    utime=calcperiodHBM(sys,U(:,ia),idisp,time);
    Uplot(:,jj,ii)=utime;
  end
end


% Coloriage de la surface
% =======================

Ucol=Uplot(:,:,3); % En fonction de la position sur l'axe des z
Ucol=ones(Nech,1)*(1:Na); % En fonction de l'indice de l'orbite périodique

% En fonction de l'abscisse curviligne ss sur la ligne (u1,u2,u3)
% définie ci-après.
ind=floor(Nech/2);
u1=Uplot(ind,:,1);
u2=Uplot(ind,:,2);
u3=Uplot(ind,:,3);
ss=sqrt((u1(2:end)-u1(1:end-1)).^2+(u2(2:end)-u2(1:end-1)).^2+(u3(2:end)-u3(1:end-1)).^2);
ss=cumsum(ss);
ss=[0 ss];
Ucol=ones(Nech,1)*ss;


% Affichage
% =========

Hand=surf(Uplot(:,:,1),Uplot(:,:,2),Uplot(:,:,3),Ucol);
set(Hand,'Edgecolor','none')

axis tight
view(-50,30)
lighting gouraud
Hchild=get(gca,'Children');

for ii=1:length(Hchild)
  hc=Hchild(ii);
  if strcmp(get(hc,'type'),'light')
    delete(hc);
  end
end
camlight('left');

grid on

xlabel(['z_',num2str(Idisp(1))])
ylabel(['z_',num2str(Idisp(2))])
zlabel(['z_',num2str(Idisp(3))])
view(-40,40)
box on


