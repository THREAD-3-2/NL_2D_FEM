gfunction [sys] = get_operators_FEM(sys)
%function [sys] = get_operators(obj)
%This function compute the operator used by Diamanlab from the quadratic
%recast of the equations.

neq_tot = sys.neq_tot;
neq     = sys.neq;
neq_aux = sys.neq_aux;
ninc    = sys.ninc;

Ndof    = sys.parameters.Ndof;
Nele    = sys.parameters.Nele;
Npgele  = sys.parameters.Npgele;
NUauxpg = sys.parameters.NUauxpg;
Npg     = sys.parameters.Npg;
ipdof   = sys.parameters.ipdof;
Nactdof = sys.parameters.Nactdof;
iactdof = sys.parameters.iactdof;

Gpg      =sys.parameters.Gpg;
detJ     =sys.parameters.detJ;
DOFELE   =sys.parameters.DOFELE;
weightG  =sys.parameters.weightG;

Fext     =sys.parameters.Fext;

equations = sys.equations;
eps = 1e-14;

% compute the tensor list of equation.m , stored in sys.iC,sys.vC,sys.iL, ...
sys=get_operator(sys);     


% for storing the lists of main equations(assumed linear)
iL_m = cell(Npg,1); jL_m = cell(Npg,1); vL_m = cell(Npg,1);

% for storing the lists of the first auxiliary equations (assumed linear)
%  gradients and values of the main variables at Gauss point
iL_g1 = cell(Npg,1); jL_g1 = cell(Npg,1); vL_g1 = cell(Npg,1);

%  for storing the lists of others auxiliary equations (at Gauss point)
iC_g = cell(Npg,1); vC_g = cell(Npg,1);
iL_g = cell(Npg,1); jL_g = cell(Npg,1); vL_g = cell(Npg,1);
iQ_g = cell(Npg,1); jQ_g = cell(Npg,1); kQ_g = cell(Npg,1); vQ_g=cell(Npg,1);

idL_g = cell(Npg,1); jdL_g = cell(Npg,1); vdL_g = cell(Npg,1);
idQ_g = cell(Npg,1); jdQ_g = cell(Npg,1); kdQ_g = cell(Npg,1); vdQ_g=cell(Npg,1);

ifunct_g= cell(Npg,1) ; jfunct_g= cell(Npg,1) ; kfunct_g= cell(Npg,1) ; funct_g= cell(Npg,1);

for ie=1:Nele  
  for ig=1:Npgele
    igg= (ie-1)*Npgele + ig;       % global gauss point number
    %  first auxiliary equations (gradients and values of main variables at
    %  Gaus points
    Rdecal = Nactdof + (igg-1)*NUauxpg;
    Udecal = Nactdof + 1 + (igg-1)*NUauxpg;
    
    iL_g1{igg} = repmat(Rdecal+1:Rdecal+4,1,Ndofele)';
    jL_g1{igg} = reshape(repmat(DOFELE(ie,:),Ngrad,1),Ndofele*Ngrad,1);
    matL_g1 = -Gpg(:,:,igg);
    vL_g1{igg} = matL_g1(:);
    
    % others auxiliary equations
   
    iC_g{igg}=sys.iC + Rdecal ; vC_g=sys.vC ;
    iL_g{igg}=sys.iL + Rdecal ; jL_g{igg}=sys.jL + Udecal ; vL_g{igg}=sys.vL ;
    iQ_g{igg}=sys.iQ + Rdecal ; jQ_g{igg}=sys.jQ + Udecal ; kQ_g{igg}=sys.kQ + Udecal ; vQ_g{igg}=sys.vQ ;
    
    idL_g{igg}=sys.idL + Rdecal ; jdL_g{igg}=sys.jdL + Udecal ; vdL_g{igg}=sys.vdL ;
    idQ_g{igg}=sys.idQ + Rdecal ; jdQ_g{igg}=sys.jdQ + Udecal ; kdQ_g{igg}=sys.kdQ + Udecal ; vL_g{igg}=sys.vdQ ;
    
    ifunct_g{igg}= sys.ifunct + Rdecal ; jfunct_g{igg}= sys.jfunct + Udecal ; kfunct_g{igg}= sys.kfunct + Udecal ; funct_g{igg}=  sys.funct ;
  end  
end  

tic
%valL5=sparse(Npg*NUauxpg);
valL5=kron(speye(Npg),valLQ8R);

%figure(11)
%spy(valL5)


valL6=sparse(Ndof,Npg*NUauxpg);
valL7=sparse(Npg*NUauxpg,Ndof);
iL6 = cell(Nele*Npgele,1);
jL6 = cell(Nele*Npgele,1);
vL6 = cell(Nele*Npgele,1);
iL7 = cell(Nele*Npgele,1);
jL7 = cell(Nele*Npgele,1);
vL7 = cell(Nele*Npgele,1);

for ie=1:Nele  
  for ig=1:Npgele
    igg= (ie-1)*Npgele + ig;       % global gauss point number
    iL6{igg} = repmat(DOFELE(ie,:),1,4)';
    jL6{igg} = reshape(repmat(NUauxpg*(igg-1)+8:NUauxpg*(igg-1)+11,16,1),16*4,1);
    matL6 = Gpg(:,:,igg)'*detJ(igg)*weightG(ig);
    vL6{igg} = matL6(:);
    
    iL7{igg} = repmat(NUauxpg*(igg-1)+1:NUauxpg*(igg-1)+4,1,16)';
    jL7{igg} = reshape(repmat(DOFELE(ie,:),4,1),16*4,1);
    matL7 = -Gpg(:,:,igg);
    vL7{igg} = matL7(:);
    
    %valL6(DOFELE(ie,:),NUauxpg*(igg-1)+8:NUauxpg*(igg-1)+11)=Gpg(:,:,igg)'*detJ(igg)*weightG(ig); 
    %valL7(NUauxpg*(igg-1)+1:NUauxpg*(igg-1)+4,DOFELE(ie,:))=-Gpg(:,:,igg);
  end  
end    

valL6 = sparse(cell2mat(iL6),cell2mat(jL6),cell2mat(vL6),Ndof,Npg*NUauxpg);
valL7 = sparse(cell2mat(iL7),cell2mat(jL7),cell2mat(vL7),Npg*NUauxpg,Ndof);

valL6(ipdof,:)=[ ];
%valL6=sparse(valL6);
valL7(:,ipdof)=[ ];
%valL7=sparse(valL7);
% figure(12)
% spy(valL6)
% figure(13)
% spy(valL7)

valL=[  sparse(neq,neq) , -Fext, valL6 ;  valL7 , sparse(neq_aux,1) , valL5 ];
%figure(12)
%spy(valL)
toc


% iQ,jQ,kQ and vQ for all the Gauss point
    tic
    Nval= size(mvQele,1);
    
    mvQ3=zeros(Nval,Npg);
    miQ3=zeros(Nval,Npg);
    mjQ3=zeros(Nval,Npg);
    mkQ3=zeros(Nval,Npg);
    
    for i=1:Npg
    mvQ3(:,i)=mvQele;
    miQ3(:,i)=miQele+(i-1)*NUauxpg+neq;
    mjQ3(:,i)=mjQele+(i-1)*NUauxpg+neq+1;
    mkQ3(:,i)=mkQele+(i-1)*NUauxpg+neq+1;
    end 
    
   mvQ4 = reshape(mvQ3,Nval*Npg,1);
   miQ4 = reshape(miQ3,Nval*Npg,1);
   mjQ4 = reshape(mjQ3,Nval*Npg,1); 
   mkQ4 = reshape(mkQ3,Nval*Npg,1); 
   toc 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% constant operator
%  Z = sparse(ninc,1);
%  valc = equations(obj,Z);
% valc(abs(valc)<eps)=0;
% [obj.iC,~,obj.vC] = find(valc);

sizeC=size(sys.vC,1)

%% linear operator
%valL(abs(valL)<eps)=0;
[sys.iL,sys.jL,sys.vL] = find(valL);

sizeL=size(sys.vL,1)

%% quadratic operator
% obj.vQ = cell2mat(vQ);
% obj.iQ = cell2mat(iQ);
% obj.jQ = cell2mat(jQ);
% obj.kQ = cell2mat(kQ);

sys.vQ = mvQ4;
sys.iQ = miQ4;
sys.jQ = mjQ4;
sys.kQ = mkQ4;

sizeQ=size(sys.vQ,1)


end
