function [obj] = get_operators_vec_FE2(obj)
%function [sys] = get_operators(obj)
%This function compute the operator used by Diamanlab from the quadratic
%recast of the equations.

neq_tot = obj.neq_tot;
neq     = obj.neq;
neq_aux = obj.neq_aux;
ninc    = obj.ninc;
Ndof    = obj.parameters.Ndof;
Nele    = obj.parameters.Nele;
Npgele  = obj.parameters.Npgele;
NUauxpg = obj.parameters.NUauxpg;
Npg     = obj.parameters.Npg;
ipdof   = obj.parameters.ipdof;
Nactdof = obj.parameters.Nactdof;
iactdof = obj.parameters.iactdof;

Gpg      =obj.parameters.Gpg;
detJ     =obj.parameters.detJ;
DOFELE   =obj.parameters.DOFELE;
weightG  =obj.parameters.weightG;

Fext     =obj.parameters.Fext;

equations = obj.equations;
eps = 1e-14;

obj.iC = []; obj.vC = [];
obj.iL = []; obj.jL = []; obj.vL = [];
obj.iQ = []; obj.jQ=[]; obj.kQ=[]; obj.vQ=[];

tic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% constant, linear and quadratic operator for a single element
% for Q8R :  Ndof = 16   NUaux = 4*11 =44  
  function [L_Qsym,dL_dQ] = eq_bilin_ele(sys,U,V,varargin)
        if nargin == 3
            [Rout1] = equations_aux_Q8R(sys,V+U);
            [Rout2] = equations_aux_Q8R(sys,V-U);
            L_Qsym = (Rout1-Rout2)/2;
        else
            [Rout1,dRout1] = equations_aux_Q8R(sys,V+U,varargin{:});
            [Rout2,dRout2] = equations_aux_Q8R(sys,V-U,varargin{:});
            L_Qsym = (Rout1-Rout2)/2;
            dL_dQ = (dRout1-dRout2)/2;
        end
 
        return
    end
    fun_eq_bilin_ele = @eq_bilin_ele;

nincele=11;
Idele=speye(nincele);
valLQ8R=zeros(nincele,nincele);
vQele = cell(nincele,1);
iQele = cell(nincele,1);
jQele = cell(nincele,1);
kQele = cell(nincele,1);

for j=1:nincele
    
    valLQ8R(:,j)=(equations_aux_Q8R(obj,Idele(:,j))-equations_aux_Q8R(obj,-Idele(:,j)))/2 ;  
    
    valQ = zeros(nincele,nincele);
        for k=1:nincele
           valQ(:,k) = feval(fun_eq_bilin_ele,obj,Idele(:,j),Idele(:,k));
        end
        valQ = (valQ - valLQ8R(:,j))/2;
        valQ(abs(valQ)<eps)=0;
        [ind_iQ,ind_kQ,ind_vQ] = find(valQ);
        vQele{j} = ind_vQ;
        iQele{j} = ind_iQ;
        jQele{j} = j*ones(size(ind_iQ));
        kQele{j} = ind_kQ;
end
%% constant operator  
valCele=equations_aux_Q8R(obj,zeros(nincele,1));
valCele(abs(valCele)<eps)=0;
%[miCele,~,mvCele] = find(valCele);

%% linear operator
valLQ8R(abs(valLQ8R)<eps)=0;
%[miLele,mjLele,mvLele] = find(valLQ8R) ;

%% quadratic operator
mvQele = cell2mat(vQele);
miQele = cell2mat(iQele);
mjQele = cell2mat(jQele);
mkQele = cell2mat(kQele);
toc


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

sizeC=size(obj.vC,1)

%% linear operator
%valL(abs(valL)<eps)=0;
[obj.iL,obj.jL,obj.vL] = find(valL);

sizeL=size(obj.vL,1)

%% quadratic operator
% obj.vQ = cell2mat(vQ);
% obj.iQ = cell2mat(iQ);
% obj.jQ = cell2mat(jQ);
% obj.kQ = cell2mat(kQ);

obj.vQ = mvQ4;
obj.iQ = miQ4;
obj.jQ = mjQ4;
obj.kQ = mkQ4;

sizeQ=size(obj.vQ,1)


end
