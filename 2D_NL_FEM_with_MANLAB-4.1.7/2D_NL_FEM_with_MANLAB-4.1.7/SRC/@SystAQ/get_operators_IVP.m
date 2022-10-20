function [sys] = get_operators_IVP(sys,equations,nR,nU,iR,iU)
% compute the index lists for the sparce tensor C,L,Q,dL,dQ,funct
%
% 
%  Let R(U) be a system of nR equations in nU variables 
%  where each equation Ri is either
%   - a quadratic expression 
%   - a functional expression like Uj= f(Uk) where f is a function
%
% 
%
%  let dR be a vector of nR quadratic differential expression where
%  each equation dRi is either
%   - empty
%   - the differential form of Uj= f(Uk) which is written under the form
%     dR :=  dUj -
%  
%  
%   let [ R, dR ]=equations(sys, U, dU) be a function that return the value
%   of R(U) and dR(U,dU) 
%
%   get_operator compute the sparse tensor C, L, Q and funct representing R(U) 
%   and  the dL dQ sparse tensor representing dR
%   
%
%
%   
eps = 1e-14;  %tolerance for erasing small terms

Znul = sparse(nU,1);
Id = speye(nU);

% working variables
valL = zeros(nR,nU);  valdL= zeros(nR,nU);

vQ = cell(nU,1); iQ = cell(nU,1); jQ = cell(nU,1); kQ = cell(nU,1);
vdQ = cell(nU,1); idQ = cell(nU,1); jdQ = cell(nU,1); kdQ = cell(nU,1);
    
    %% get the index of function using  dR  
    [~,dRtest] = equations(sys,randn(nU,1),randn(nU,1));
    if strcmp(sys.writing,'standard') 
     % in the standard case, diff equa are allowed in R, but they 
     % are not associated to a function, function are allowed only in Raux.
     dRtest(1:sys.neq)=0;  % remove diff equa detection in R 
    end
    ind_diff = find(dRtest);
    %if numel(dRtest) == neq_aux; ind_diff = ind_diff+neq; end    
    
    %% constant operator C
    valc = equations(sys,Znul,Znul);
    valc(ind_diff) = 0;
    valc(abs(valc)<eps)=0;
    [iC,~,vC] = find(valc);

    %% linear and quadratic operators L,Q,dL,dQ
      
    %if numel(dRtest) == neq_aux; valdL = zeros(neq_aux,ninc);else; valdL = zeros(neq_tot,ninc); end
    for j=1:nU
        
        valL(:,j) = (equations(sys,Id(:,j),Znul)-equations(sys,-Id(:,j),Znul))/2;
        [~,valdL(:,j)] = equations(sys,Znul,Id(:,j));
        
        valQ = zeros(nR,nU);
        %if numel(dRtest) == neq_aux; valdQ = zeros(neq_aux,ninc);else; valdQ = zeros(neq_tot,ninc); end
        valdQ = zeros(nR,nU);
        for k=1:nU
            valQ(:,k) = (equations(sys,Id(:,k)+Id(:,j),Znul)-equations(sys,Id(:,k)-Id(:,j),Znul))/2;
            [~,valdQ(:,k)] = equations(sys,Id(:,k),Id(:,j));
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
    
    valL(ind_diff,:) = 0;
    valL(abs(valL)<eps)=0;
    [iL,jL,vL] = find(valL);
    valdL(abs(valdL)<eps)=0;
    [idL,jdL,vdL] = find(valdL);
    
    
    iQ = cell2mat(iQ(~cellfun('isempty',iQ))); jQ = cell2mat(jQ(~cellfun('isempty',jQ)));
    kQ = cell2mat(kQ(~cellfun('isempty',kQ))); vQ = cell2mat(vQ(~cellfun('isempty',vQ)));
    
    idQ= cell2mat(idQ(~cellfun('isempty',idQ))); jdQ= cell2mat(jdQ(~cellfun('isempty',jdQ)));
    kdQ= cell2mat(kdQ(~cellfun('isempty',kdQ))); vdQ= cell2mat(vdQ(~cellfun('isempty',vdQ)));
    
    %if numel(dRtest) == neq_aux; sys.idQ = sys.idQ+neq; sys.idL = sys.idL+neq; end

    %% get the indexes for each function  cos, sin , exp, ... in the R vector
    
    [Rtest,~] = equations(sys,ones(nU,1)/2,Znul);
    
    
    ifunct =ind_diff; % equation lines  with : U(k) = funct( U(j))
    ndiff=size(ind_diff,1);
    jfunct =zeros(ndiff,1); kfunct =zeros(ndiff,1); funct =cell(ndiff,1);
    % retrieves which function it is by evaluating f(1/2),  
    % then retrieves j, k using dQ, dL
    for i=1:ndiff
       
        if abs(Rtest(ifunct(i))-(1/2-exp(1/2)))<eps; funct{i}=@exp; end
        if abs(Rtest(ifunct(i))-(1/2-log(1/2)))<eps; funct{i}=@log; end
        
        if abs(Rtest(ifunct(i))-(1/2-cos(1/2)))<eps; funct{i}=@cos; end
        if abs(Rtest(ifunct(i))-(1/2-sin(1/2)))<eps; funct{i}=@sin; end
        if abs(Rtest(ifunct(i))-(1/2-tan(1/2)))<eps; funct{i}=@tan; end
        if abs(Rtest(ifunct(i))-(1/2-acos(1/2)))<eps; funct{i}=@acos; end
        if abs(Rtest(ifunct(i))-(1/2-asin(1/2)))<eps; funct{i}=@asin; end
        if abs(Rtest(ifunct(i))-(1/2-atan(1/2)))<eps; funct{i}=@atan; end
        if abs(Rtest(ifunct(i))-(1/2-cosh(1/2)))<eps; funct{i}=@cosh; end
        if abs(Rtest(ifunct(i))-(1/2-sinh(1/2)))<eps; funct{i}=@sinh; end
        if abs(Rtest(ifunct(i))-(1/2-tanh(1/2)))<eps; funct{i}=@tanh; end
       
        jfunct(i)=kdQ(find(idQ==ifunct(i)));    
        kfunct(i)=jdL(find(idL==ifunct(i)));
    end  
    
% Actualize sys lists
    
switch sys.writing
    case 'standard'      % standard case  full system in @equatiuons
     sys.iC  = iC ; sys.vC = vC ;
     sys.iL  = iL ; sys.jL = jL ; sys.vL = vL ;
     sys.iQ  = iQ ; sys.jQ = jQ ; sys.kQ  =kQ ; sys.vQ = vQ;
     
     sys.idL = idL ; sys.jdL = jdL ; sys.vdL= vdL;
     sys.idQ = idQ ; sys.jdQ = jdQ ; sys.kdQ = kdQ ;sys.vdQ= vdQ;
     
     sys.ifunct = ifunct ; sys.jfunct = jfunct ; sys.kfunct = kfunct ; sys.funct = funct;
     
    case 'subsystems'     % subsystem case : add  @equations  to sys lists
     sys.iC=[sys.iC;iR(iC)]; sys.vC=[sys.vC;vC];
     sys.iL=[sys.iL;iR(iL)]; sys.jL=[sys.jL;iU(jL)]; sys.vL=[sys.vL;vL];
     sys.iQ=[sys.iQ;iR(iQ)]; sys.jQ=[sys.jQ;iU(jQ)]; sys.kQ=[sys.kQ;iU(kQ)]; sys.vQ=[sys.vQ;vQ]; 
     
     sys.idL=[sys.idL;iR(idL)]; sys.jdL=[sys.jdL;iU(jdL)]; sys.vdL=[sys.vdL;vdL];
     sys.idQ=[sys.idQ;iR(idQ)]; sys.jdQ=[sys.jdQ;iU(jdQ)]; sys.kdQ=[ sys.kdQ;iU(kdQ)]; sys.vdQ=[sys.vdQ;vdQ];
     
     sys.ifunct=[sys.ifunct;iR(ifunct)]; sys.jfunct=[sys.jfunct;iU(jfunct)]; sys.kfunct=[sys.kfunct;iU(kfunct)];
     sys.funct=[sys.funct;funct] ;
end 

end





