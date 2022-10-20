function [sys] = get_operators(sys)
%function [sys] = get_operators(sys)
%This function compute the operator used by Manlab from the quadratic
%recast of the equations using polarization formulaes.

nz_tot = sys.nz_tot;
nz_aux = sys.nz_aux;
nz = sys.nz;
equations = sys.equations;
eps = 1e-14;

    function [L_Qsym,Fsym] = eq_bilin(U,V,varargin)
        if nargout == 1
            [Rout1] = equations(cell_argin_eq{:},V+U,varargin{:});
            [Rout2] = equations(cell_argin_eq{:},V-U,varargin{:});
            % Polarization formula
            L_Qsym = (Rout1-Rout2)/2;
        else
            [Rout1,~,Forced1] = equations(cell_argin_eq{:},V+U,varargin{:});
            [Rout2,~,Forced2] = equations(cell_argin_eq{:},V-U,varargin{:});
            % Polarization formula
            L_Qsym = (Rout1-Rout2)/2;
            Fsym = (Forced1-Forced2)/2;
        end
        
        return
    end

sys.ic0 = []; sys.vc0 = [];
sys.ic1 = []; sys.vc1 = [];
sys.ic2 = []; sys.vc2 = [];

sys.il0 = []; sys.jl0 = []; sys.vl0 = [];
sys.il1 = []; sys.jl1 = []; sys.vl1 = [];

sys.id  = []; sys.jd  = []; sys.vd  = [];
sys.idd  = []; sys.jdd  = []; sys.vdd  = [];
sys.id1  = []; sys.jd1  = []; sys.vd1  = [];
sys.idl = []; sys.jdl = []; sys.vdl = [];

sys.iq =[]; sys.jq =[]; sys.kq =[]; sys.vq =[];
sys.idq =[]; sys.jdq =[]; sys.kdq =[]; sys.vdq =[];

sys.iforce0 = []; sys.hforce0 = []; sys.vforce0 = [];
sys.iforce1 = []; sys.hforce1 = []; sys.vforce1 = [];
sys.iforce2 = []; sys.hforce2 = []; sys.vforce2 = [];

H = sys.H;
Z = sparse(nz_tot+1,1);
Id = speye(nz_tot+1);

fun_eq_bilin = @eq_bilin;

vq = cell(nz_tot,1);
iq = cell(nz_tot,1);
jq = cell(nz_tot,1);
kq = cell(nz_tot,1);
vdq = cell(nz_tot,1);
idq = cell(nz_tot,1);
jdq = cell(nz_tot,1);
kdq = cell(nz_tot,1);

cell_argin_eq = cell(2,1);
cell_argin_eq{1} = sys;
cell_argin_eq{2} = 0;

[~,dRtest] = feval(equations,cell_argin_eq{:},randn(nz_tot+1,1),randn(nz_tot,1),Z);
ind_diff = find(dRtest)+nz;

%% constant operators C0 and C1
valc0 = feval(equations,cell_argin_eq{:},Z,Z,Z);
valc1 = feval(fun_eq_bilin,Id(:,end),Z,Z,Z);
valc2 = (feval(fun_eq_bilin,Id(:,end),Id(:,end),Z,Z) - valc1)/2;
valc0(ind_diff) = 0;
valc1(ind_diff) = 0;
valc2(ind_diff) = 0;
valc0(abs(valc0)<eps)=0;
valc1(abs(valc1)<eps)=0;
valc2(abs(valc2)<eps)=0;

%% linear/derivation operators L0, L1, d1, dd
val0 = zeros(nz_tot,nz_tot);
val1 = zeros(nz_tot,nz_tot);
vald = zeros(nz_tot,nz_tot);
valdl = zeros(nz_aux,nz_tot);
valdd = zeros(nz_tot,nz_tot);
vald1 = zeros(nz_tot,nz_tot);
parfor j=1:nz_tot
    val0(:,j) = feval(fun_eq_bilin,Id(:,j),Z,Z,Z);
    val1(:,j) = feval(fun_eq_bilin,Id(:,j),Id(:,end),Z,Z);
    valdd(:,j) = feval(equations,cell_argin_eq{:},Z,Z,Id(:,j));
    [vald(:,j),valdl(:,j)] = feval(equations,cell_argin_eq{:},Z,Id(:,j),Z);
    vald1(:,j) = feval(equations,cell_argin_eq{:},Id(:,end),Id(:,j),Z);
end

val0(ind_diff,:) = 0;
val1 = val1 - val0;
val1(ind_diff,:) = 0;
vald = vald - valc0;
vald(ind_diff,:) = 0;
vald1 = vald1 - valc0 - valc1 - valc2 - vald;
vald1(ind_diff,:) = 0;
valdd = valdd - valc0;
valdd(ind_diff,:) = 0;

val0(abs(val0)<eps)=0;
val1(abs(val1)<eps)=0;
vald(abs(vald)<eps)=0;
valdl(abs(valdl)<eps)=0;
valdd(abs(valdd)<eps)=0;
vald1(abs(vald1)<eps)=0;

%% Quadratic operator Q
parfor j=1:nz_tot
    valq = zeros(nz_tot,nz_tot);
    valdq = zeros(nz_aux,nz_tot);
    for k=1:nz_tot
        valq(:,k) = feval(fun_eq_bilin,Id(:,j),Id(:,k),Z,Z);
        [~,valdq(:,k)] = feval(equations,cell_argin_eq{:},Id(:,j),Id(:,k),Z);
    end
    valq = (valq - val0(:,j))/2;
    valq(ind_diff,:) = 0;
    valq(abs(valq)<eps)=0;
    [ind_iq,ind_kq,ind_vq] = find(valq);
    vq{j} = ind_vq;
    iq{j} = ind_iq;
    jq{j} = j*ones(size(ind_iq));
    kq{j} = ind_kq;
    
    valdq = valdq - valdl;
    valdq(abs(valdq)<eps)=0;
    [ind_idq,ind_kdq,ind_vdq] = find(valdq);
    vdq{j} = ind_vdq;
    idq{j} = ind_idq;
    jdq{j} = j*ones(size(ind_idq));
    kdq{j} = ind_kdq;
end

%% Forced terms
t = 0:2*pi/(2*H+1):2*pi*(2*H/(2*H+1));
cell_argin_eq{2} = t;
[~,~,forced_terms0] = feval(equations,sys,t,Z,Z,Z);
[~,forced_terms1] = feval(fun_eq_bilin,Id(:,end),Z,Z,Z);
[~,forced_terms2] = feval(fun_eq_bilin,Id(:,end),Id(:,end),Z,Z);

FFTforce0 = fft(forced_terms0)/(2*H+1);
FFTforce1 = fft(forced_terms1)/(2*H+1);
FFTforce2 = fft( (forced_terms2-forced_terms1)/2 )/(2*H+1);

valforce0 = sys.compl_to_real(FFTforce0(1:H+1,:));
valforce1 = sys.compl_to_real(FFTforce1(1:H+1,:));
valforce2 = sys.compl_to_real(FFTforce2(1:H+1,:));

valforce0(abs(valforce0)<eps)=0;
valforce1(abs(valforce1)<eps)=0;
valforce2(abs(valforce2)<eps)=0;

% Forced terms
[sys.iforce0,sys.hforce0,sys.vforce0] = find(valforce0.');
[sys.iforce1,sys.hforce1,sys.vforce1] = find(valforce1.');
[sys.iforce2,sys.hforce2,sys.vforce2] = find(valforce2.');

% Constant operators
[sys.ic0,~,sys.vc0] = find(valc0);
[sys.ic1,~,sys.vc1] = find(valc1);
[sys.ic2,~,sys.vc2] = find(valc2);

% Linear operators
[sys.il0,sys.jl0,sys.vl0] = find(val0);
[sys.il1,sys.jl1,sys.vl1] = find(val1);
[sys.idl,sys.jdl,sys.vdl] = find(valdl);
sys.idl = sys.idl+nz;

% Derivation operators
[sys.id,sys.jd,sys.vd] = find(vald);
[sys.idd,sys.jdd,sys.vdd] = find(valdd);
[sys.id1,sys.jd1,sys.vd1] = find(vald1);

% Quadratic operators
sys.vq = cell2mat(vq);
sys.iq = cell2mat(iq);
sys.jq = cell2mat(jq);
sys.kq = cell2mat(kq);

sys.vdq = cell2mat(vdq);
sys.idq = cell2mat(idq)+nz;
sys.jdq = cell2mat(jdq);
sys.kdq = cell2mat(kdq);

end