function [sys] = get_operators(sys)
%function [sys] = get_operators(obj)
%This function compute the operator used by Diamanlab from the quadratic
%recast of the equations.

neq_tot = sys.neq_tot;
neq_aux = sys.neq_aux;
neq = sys.neq;
ninc = sys.ninc;
equations = sys.equations;
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

sys.iC = []; sys.vC = [];
sys.iL = []; sys.jL = []; sys.vL = [];
sys.iQ = []; sys.jQ=[]; sys.kQ=[]; sys.vQ=[];

Z = sparse(ninc,1);
Id = speye(ninc);

valL = zeros(neq_tot,ninc);
fun_eq_bilin = @eq_bilin;

vQ = cell(ninc,1);
iQ = cell(ninc,1);
jQ = cell(ninc,1);
kQ = cell(ninc,1);

if nargin(equations) == 2
    
    %% constant operator
    valc = equations(sys,Z);
    
    %% linear and quadratic operators
    
    %%% Parallelization of the construction for more efficiency.
    for j=1:ninc
        
        valL(:,j) = feval(fun_eq_bilin,sys,Id(:,j),Z);
        
        valQ = zeros(neq_tot,ninc)
        for k=1:ninc
            valQ(:,k) = feval(fun_eq_bilin,sys,Id(:,j),Id(:,k))
        end
        valQ = (valQ - valL(:,j))/2;
        valQ(abs(valQ)<eps)=0;
        [ind_iQ,ind_kQ,ind_vQ] = find(valQ);
        vQ{j} = ind_vQ;
        iQ{j} = ind_iQ;
        jQ{j} = j*ones(size(ind_iQ));
        kQ{j} = ind_kQ;
    end
    
elseif nargin(equations) == 3
    
    [~,dRtest] = equations(sys,randn(ninc,1),randn(ninc,1));
    ind_diff = find(dRtest);
    if numel(dRtest) == neq_aux; ind_diff = ind_diff+neq; end

    sys.idL = []; sys.jdL = []; sys.vdL = [];
    sys.idQ = []; sys.jdQ=[]; sys.kdQ=[]; sys.vdQ=[];
    
    %% constant operator
    valc = equations(sys,Z,Z);
    valc(ind_diff) = 0;

    %% linear and quadratic operators
    vdQ = cell(ninc,1);
    idQ = cell(ninc,1);
    jdQ = cell(ninc,1);
    kdQ = cell(ninc,1);
    
    if numel(dRtest) == neq_aux; valdL = zeros(neq_aux,ninc);
    else; valdL = zeros(neq_tot,ninc); end
    
    %%% Parallelization of the construction for more efficiency.
    for j=1:ninc
        
        valL(:,j) = feval(fun_eq_bilin,sys,Id(:,j),Z,Z);
        [~,valdL(:,j)] = feval(equations,sys,Z,Id(:,j));
        
        valQ = zeros(neq_tot,ninc);
        if numel(dRtest) == neq_aux; valdQ = zeros(neq_aux,ninc);
        else; valdQ = zeros(neq_tot,ninc); end
        for k=1:ninc
            valQ(:,k) = feval(fun_eq_bilin,sys,Id(:,j),Id(:,k),Z);
            [~,valdQ(:,k)] = feval(equations,sys,Id(:,k),Id(:,j));
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
    [sys.idL,sys.jdL,sys.vdL] = find(valdL);
    
    sys.vdQ = cell2mat(vdQ(~cellfun('isempty',vdQ)));
    sys.idQ = cell2mat(idQ(~cellfun('isempty',idQ)));
    sys.jdQ = cell2mat(jdQ(~cellfun('isempty',jdQ)));
    sys.kdQ = cell2mat(kdQ(~cellfun('isempty',kdQ)));
    if numel(dRtest) == neq_aux; sys.idQ = sys.idQ+neq; sys.idL = sys.idL+neq; end

    valL(ind_diff,:) = 0;
end

%% constant operator
valc(abs(valc)<eps)=0;
[sys.iC,~,sys.vC] = find(valc);

%% linear operator
valL(abs(valL)<eps)=0;
[sys.iL,sys.jL,sys.vL] = find(valL);

%% quadratic operator
sys.vQ = cell2mat(vQ(~cellfun('isempty',vQ)));
sys.iQ = cell2mat(iQ(~cellfun('isempty',iQ)));
sys.jQ = cell2mat(jQ(~cellfun('isempty',jQ)));
sys.kQ = cell2mat(kQ(~cellfun('isempty',kQ)));
end