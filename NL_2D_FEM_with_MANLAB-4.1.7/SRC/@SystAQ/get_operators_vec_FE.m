function [obj] = get_operators_vec_FE(obj)
%function [sys] = get_operators(obj)
%This function compute the operator used by Diamanlab from the quadratic
%recast of the equations.

neq_tot = obj.neq_tot;
neq     = obj.neq;
neq_aux = obj.neq_aux;
ninc    = obj.ninc;
NUauxpg = obj.parameters.NUauxpg;
Npg     = obj.parameters.Npg;

equations = obj.equations;
eps = 1e-14;

    function [L_Qsym,dL_dQ] = eq_bilin(sys,U,V,varargin)
        if nargin == 3
            [Rout1] = equations(sys,V+U);
            [Rout2] = equations(sys,V-U);
            L_Qsym = (Rout1-Rout2)/2;
        else
            [Rout1,dRout1] = equations(sys,V+U,varargin{:});
            [Rout2,dRout2] = equations(sys,V-U,varargin{:});
            L_Qsym = (Rout1-Rout2)/2;
            dL_dQ = (dRout1-dRout2)/2;
        end
 
        return
    end

obj.iC = []; obj.vC = [];
obj.iL = []; obj.jL = []; obj.vL = [];
obj.iQ = []; obj.jQ=[]; obj.kQ=[]; obj.vQ=[];

Z = sparse(ninc,1);
Id = speye(ninc);

valL = sparse(neq_tot,ninc);
valL2 = sparse(neq_tot,ninc);
fun_eq_bilin = @eq_bilin;

vQ = cell(ninc,1);
iQ = cell(ninc,1);
jQ = cell(ninc,1);
kQ = cell(ninc,1);

vQ2 = cell(ninc,1);
iQ2 = cell(ninc,1);
jQ2 = cell(ninc,1);
kQ2 = cell(ninc,1);


if nargin(equations) == 2
    
    %% constant operator
    valc = equations(obj,Z);
    
    %% linear and quadratic operators
    
    %%% Parallelization of the construction for more efficiency.
%     for j=1:ninc
%         
%         valL(:,j) = feval(fun_eq_bilin,obj,Id(:,j),Z);
%         
%         valQ = zeros(neq_tot,ninc);
%         for k=1:ninc
%             valQ(:,k) = feval(fun_eq_bilin,obj,Id(:,j),Id(:,k));
%         end
%         valQ = (valQ - valL(:,j))/2;
%         valQ(abs(valQ)<eps)=0;
%         [ind_iQ,ind_kQ,ind_vQ] = find(valQ);
%         vQ{j} = ind_vQ;
%         iQ{j} = ind_iQ;
%         jQ{j} = j*ones(size(ind_iQ));
%         kQ{j} = ind_kQ;
%     end
%     % sparce tensor list
%     mvQ = cell2mat(vQ)
%     miQ = cell2mat(iQ)
%     mjQ = cell2mat(jQ)
%     mkQ = cell2mat(kQ) 
   
   
    %  if R(U)=C + L(U) + Q(U,U)
    %  then R(U)-R(-U)= 2*L(U)
    %  linear term captur
   
%    for j=1:ninc
%        valL2(:,j)=(equations(obj,Id(:,j))-equations(obj,-Id(:,j)))/2 ;
%    end
     %tic
     %valL=(equations(obj,Id)-equations(obj,-Id))/2 ;
     %toc
     
     tic
     valL2(:,1:neq+1)=(equations(obj,Id(:,1:neq+1))-equations(obj,-Id(:,1:neq+1)))/2 ;
     toc
     tic
     for i=1:NUauxpg
        ideb=neq+2+Npg*(i-1);
        ifin=neq+1+Npg*i;
        valL2(:,ideb:ifin)=(equations(obj,Id(:,ideb:ifin))-equations(obj,-Id(:,ideb:ifin)))/2 ;   
     end    
    toc
    
    spy(valL2)
    
    % iQ,jQ,kQ and vQ for the first Gauss point
    tic
    for j=1:NUauxpg
        index_j=neq+1+Npg*(j-1)+1;
        Utotj=zeros(ninc,1); Utotj(index_j)=1;
        
        valQ = zeros(neq_tot,NUauxpg);
        for k=1:NUauxpg
             index_k=neq+1+Npg*(k-1)+1;
             Utotk=zeros(ninc,1); Utotk(index_k)=1;
             valQ(:,k) = feval(fun_eq_bilin,obj,Utotj,Utotk);
        end
        valQ = (valQ - valL2(:,index_j))/2;
        valQ(abs(valQ)<eps)=0;
        [ind_iQ,ind_kQ,ind_vQ] = find(valQ);
        vQ2{j} = ind_vQ;
        iQ2{j} = ind_iQ;
        jQ2{j} = index_j*ones(size(ind_iQ));
        kQ2{j} = (ind_kQ-1)*Npg+neq+2;
    end
    mvQ2 = cell2mat(vQ2);
    miQ2 = cell2mat(iQ2);
    mjQ2 = cell2mat(jQ2);
    mkQ2 = cell2mat(kQ2);
    
    toc
    % iQ,jQ,kQ and vQ for all the Gauss point
    tic
    Nval= size(mvQ2,1);
    
    mvQ3=zeros(Nval,Npg);
    miQ3=zeros(Nval,Npg);
    mjQ3=zeros(Nval,Npg);
    mkQ3=zeros(Nval,Npg);
    
    for i=1:Npg
    mvQ3(:,i)=mvQ2;
    miQ3(:,i)=miQ2+i-1;
    mjQ3(:,i)=mjQ2+i-1;
    mkQ3(:,i)=mkQ2+i-1;
    end 
    
   mvQ4 = reshape(mvQ3,Nval*Npg,1);
   miQ4 = reshape(miQ3,Nval*Npg,1);
   mjQ4 = reshape(mjQ3,Nval*Npg,1); 
   mkQ4 = reshape(mkQ3,Nval*Npg,1); 
   toc 
    
    
    
    
elseif nargin(equations) == 3
    
    [~,dRtest] = equations(obj,randn(ninc,1),randn(ninc,1));
    ind_diff = find(dRtest) + neq;

    obj.idL = []; obj.jdL = []; obj.vdL = [];
    obj.idQ = []; obj.jdQ=[]; obj.kdQ=[]; obj.vdQ=[];
    
    %% constant operator
    valc = equations(obj,Z,Z);
        valc(ind_diff) = 0;

    %% linear and quadratic operators
    vdQ = cell(ninc,1);
    idQ = cell(ninc,1);
    jdQ = cell(ninc,1);
    kdQ = cell(ninc,1);
    
    valdL = zeros(neq_aux,ninc);
    %%% Parallelization of the construction for more efficiency.
    for j=1:ninc
        
        valL(:,j) = feval(fun_eq_bilin,obj,Id(:,j),Z,Z);
        [~,valdL(:,j)] = feval(equations,obj,Z,Id(:,j));
        
        valQ = zeros(neq_tot,ninc);
        valdQ = zeros(neq_aux,ninc);
        for k=1:ninc
            valQ(:,k) = feval(fun_eq_bilin,obj,Id(:,j),Id(:,k),Z);
            [~,valdQ(:,k)] = feval(equations,obj,Id(:,k),Id(:,j));
        end
        valQ = (valQ - valL(:,j))/2;
        valQ(ind_diff,:) = 0;
        valQ(abs(valQ)<eps)=0;
        [ind_iQ,ind_kQ,ind_vQ] = find(valQ);
        vQ{j} = ind_vQ;
        iQ{j} = ind_iQ;
        jQ{j} = j*ones(size(ind_iQ));
        kQ{j} = ind_kQ;
        
        valdQ = (valdQ - valdL(:,j));
        valdQ(abs(valdQ)<eps)=0;
        [ind_idQ,ind_kdQ,ind_vdQ] = find(valdQ);
        vdQ{j} = ind_vdQ;
        idQ{j} = ind_idQ;
        kdQ{j} = j*ones(size(ind_idQ));
        jdQ{j} = ind_kdQ;
    end
    
    valdL(abs(valdL)<eps)=0;
    [obj.idL,obj.jdL,obj.vdL] = find(valdL);
    
    obj.vdQ = cell2mat(vdQ);
    obj.idQ = cell2mat(idQ);
    obj.jdQ = cell2mat(jdQ);
    obj.kdQ = cell2mat(kdQ);
    
    % Indices augmented by neq : only Raux can be affected by differentiation.
    obj.idL = obj.idL + obj.neq;
    obj.idQ = obj.idQ + obj.neq;
    
    valL(ind_diff,:) = 0;
end


%% constant operator
valc(abs(valc)<eps)=0;
[obj.iC,~,obj.vC] = find(valc);

%% linear operator
tic
valL(abs(valL)<eps)=0;
%[obj.iL,obj.jL,obj.vL] = find(valL);
[obj.iL,obj.jL,obj.vL] = find(valL2);
toc
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
