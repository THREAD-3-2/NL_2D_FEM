function Jacobian = Jacobian_QPHBM(sys,Utot)
% Compute the Jacobian of R with respect to U

H    = sys.H;
nz_tot   = sys.nz_tot;

coef_per_var = (H+1)*4*H+1;

% Extract informations from Utot :
[Ztot,omega1,omega2,lambda,omega1sq,omega2sq,omega1omega2,lambdaomega1,lambdaomega2] = get_Ztot(sys,Utot);

Z_t0 = sum(Ztot(1:(coef_per_var-1)/2+1,:),1);

% definition of the small matrices : I and d
%
% Identity matrix coef_per_var sparse
[iI,jI,vI]=find(speye(coef_per_var)) ;

% derivation indices
cos_ind = 2:(coef_per_var-1)/2+1;
sin_ind = (coef_per_var-1)/2+2:coef_per_var;

Domega1vec = [zeros(H,1) ; reshape((1:H).*ones(2*H+1,1),2*H*H+H,1)];
Domega2vec = [(1:H)' ; reshape((-H:1:H)'.*ones(1,H),2*H*H+H,1)];
Dvec = omega1*Domega1vec+omega2*Domega2vec;
D1vec = lambdaomega1*Domega1vec+lambdaomega2*Domega2vec;
%DLvec = Domega1vec+Domega2vec;
DDvec = omega1sq*(Domega1vec.*Domega1vec) + omega2sq*(Domega2vec.*Domega2vec) + (2*omega1omega2)*(Domega1vec.*Domega2vec);

% d sparse matrix
idmat=[ 1 ,   cos_ind ,  sin_ind ]';
jdmat=[ 1 ,   sin_ind ,  cos_ind ]';
vdmat= [0;Dvec ;-Dvec ];
vd1mat=[0;D1vec;-D1vec];
vdlmat=[0;Dvec;-Dvec];

% dd sparse matrix
iddmat=idmat;
jddmat=iddmat;
vddmat=[0;-DDvec;-DDvec];

% ib and jb for b sparce matrix
[ib,jb,~]=find( sparse( ones(coef_per_var) ));

%   dRdZ =  L0 + lambda L1  + 2*Q(Z,.)
%   construct the list for the large matrix L0, L1, dl, 2*Q(Z,.) and dq(Z,.)
nl0=numel(sys.il0); iL0=zeros(coef_per_var,nl0); jL0=zeros(coef_per_var,nl0); vL0=zeros(coef_per_var,nl0);
for i=1:nl0
    iL0(:,i)=iI + (sys.il0(i)-1)*coef_per_var ;
    jL0(:,i)=jI + (sys.jl0(i)-1)*coef_per_var ;
    vL0(:,i)=vI*sys.vl0(i)  ;
end

nl1=numel(sys.il1); iL1=zeros(coef_per_var,nl1); jL1=zeros(coef_per_var,nl1); vL1=zeros(coef_per_var,nl1);
for i=1:nl1
    iL1(:,i)=iI + (sys.il1(i)-1)*coef_per_var ;
    jL1(:,i)=jI + (sys.jl1(i)-1)*coef_per_var ;
    vL1(:,i)=vI*(sys.vl1(i)*lambda) ;
end

nq=numel(sys.iq); iQL=zeros(coef_per_var^2,nq); jQL=zeros(coef_per_var^2,nq); vQL=zeros(coef_per_var^2,nq);
for i=1:nq
    % Avoid to compute the same B several times / it uses that obj.jq is non-decreasing.
    if ~(i>1 && sys.jq(i-1) == sys.jq(i))
        vb= reshape(sys.B(Ztot(:,sys.jq(i))) ,coef_per_var^2,1) ;
    end
    iQL(:,i)=ib + (sys.iq(i)-1)*coef_per_var ;
    jQL(:,i)=jb + (sys.kq(i)-1)*coef_per_var ;
    vQL(:,i)=2*vb*sys.vq(i) ;
end

nld=numel(sys.id); iD=zeros(coef_per_var,nld); jD=zeros(coef_per_var,nld); vD=zeros(coef_per_var,nld);
for i=1:nld
    iD(:,i)=idmat + (sys.id(i)-1)*coef_per_var ;
    jD(:,i)=jdmat + (sys.jd(i)-1)*coef_per_var ;
    vD(:,i)=vdmat*sys.vd(i) ;
end

nld1=numel(sys.id1); id1=zeros(coef_per_var,nld1); jd1=zeros(coef_per_var,nld1); vd1=zeros(coef_per_var,nld1);
for i=1:nld1
    id1(:,i)=idmat + (sys.id1(i)-1)*coef_per_var ;
    jd1(:,i)=jdmat + (sys.jd1(i)-1)*coef_per_var ;
    vd1(:,i)=vd1mat*sys.vd1(i) ;
end

nldd=numel(sys.idd); idd=zeros(coef_per_var,nldd); jdd=zeros(coef_per_var,nldd); vdd=zeros(coef_per_var,nldd);
for i=1:nldd
    idd(:,i)=iddmat + (sys.idd(i)-1)*coef_per_var ;
    jdd(:,i)=jddmat + (sys.jdd(i)-1)*coef_per_var ;
    vdd(:,i)=vddmat*sys.vdd(i) ;
end



% Assumption : the variable after kdq(i) is the derivative of the variable kdq(i)

nldl=numel(sys.idl); idl=zeros(coef_per_var,nldl); jdl=zeros(coef_per_var,nldl); vdl=zeros(coef_per_var,nldl);
for i=1:nldl
    idl(:,i)=idmat + (sys.idl(i)-1)*coef_per_var ;
    jdl(:,i)=jdmat + (sys.jdl(i)-1)*coef_per_var ;
    vdl(:,i)=vdlmat*sys.vdl(i) ;
end

% Assumption : the variable after kdq(i) is the derivative of the variable kdq(i)
% Assumption : the variable after kdq(i) is the derivative of the variable kdq(i)

ndq=numel(sys.idq); idq1=zeros(coef_per_var^2,ndq); jdq1=zeros(coef_per_var^2,ndq); vdq1=zeros(coef_per_var^2,ndq);
for i=1:ndq
    vb= reshape(sys.B(Ztot(:,sys.jdq(i))),coef_per_var^2,1) ;
    idq1(:,i)=ib + (sys.idq(i)-1)*coef_per_var ;
    % Assumption : the variable after kdq(i) is the derivative of the variable kdq(i)
    jdq1(:,i)=jb + (sys.kdq(i)+1-1)*coef_per_var ;
    vdq1(:,i)=vb*sys.vdq(i) ;
end

idq2=zeros(coef_per_var^2,ndq); jdq2=zeros(coef_per_var^2,ndq); vdq2=zeros(coef_per_var^2,ndq);
for i=1:ndq
    % Assumption : the variable after kdq(i) is the derivative of the variable kdq(i)
    vb= reshape(sys.B(Ztot(:,sys.kdq(i)+1)),coef_per_var^2,1) ;
    idq2(:,i)=ib + (sys.idq(i)-1)*coef_per_var ;
    jdq2(:,i)=jb + (sys.jdq(i)-1)*coef_per_var ;
    vdq2(:,i)=vb*sys.vdq(i) ;
end



idRdZ=[ reshape(iL0,nl0*coef_per_var,1) ; reshape(iL1,nl1*coef_per_var,1) ; ...
    reshape(iQL,nq*coef_per_var^2,1) ; reshape(iD,nld*coef_per_var,1); ...
    reshape(id1,nld1*coef_per_var,1) ; reshape(idd,nldd*coef_per_var,1); ...
    reshape(idq2,ndq*coef_per_var^2,1) ; reshape(idq1,ndq*coef_per_var^2,1) ; reshape(idl,nldl*(coef_per_var),1)];

jdRdZ=[ reshape(jL0,nl0*coef_per_var,1) ; reshape(jL1,nl1*coef_per_var,1) ; ...
    reshape(jQL,nq*coef_per_var^2,1) ; reshape(jD,nld*coef_per_var,1); ...
    reshape(jd1,nld1*coef_per_var,1) ; reshape(jdd,nldd*coef_per_var,1); ...
    reshape(jdq2,ndq*coef_per_var^2,1) ; reshape(jdq1,ndq*coef_per_var^2,1) ; reshape(jdl,nldl*(coef_per_var),1)];

vdRdZ=[ reshape(vL0,nl0*coef_per_var,1) ; reshape(vL1,nl1*coef_per_var,1) ; ...
    reshape(vQL,nq*coef_per_var^2,1) ; reshape(vD,nld*coef_per_var,1); ...
    reshape(vd1,nld1*coef_per_var,1) ; reshape(vdd,nldd*coef_per_var,1); ...
    reshape(vdq2,ndq*coef_per_var^2,1) ; reshape(vdq1,ndq*coef_per_var^2,1) ; reshape(vdl,nldl*(coef_per_var),1)];

dRdZ=sparse(idRdZ,jdRdZ,vdRdZ,coef_per_var*nz_tot,coef_per_var*nz_tot) ;

%   Modification for dl/dq : initial condition replace the zero-th harmonic :
dRdZ(coef_per_var*(sys.idl-1)+1,:) = 0;
for i=1:numel(sys.idl)
    dRdZ(coef_per_var*(sys.idl(i)-1)+1,coef_per_var*(sys.jdl(i)-1)+(1:(coef_per_var+1)/2)) = sys.vdl(i);
end
for i=1:numel(sys.idq)
    dRdZ(coef_per_var*(sys.idq(i)-1)+1,coef_per_var*(sys.kdq(i)-1)+(1:(coef_per_var+1)/2)) = sys.vdq(i)*Z_t0(sys.jdq(i))  +  dRdZ(coef_per_var*(sys.idq(i)-1)+1,coef_per_var*(sys.kdq(i)-1)+(1:(coef_per_var+1)/2));
end

%   dRd1ambda  = C1 + L1(Z) + 2*lambda*C2
dRdlambda =zeros(coef_per_var,nz_tot);
for i=1:numel(sys.iforce1)
    dRdlambda( (sys.hforce1(i)-1)*(2*H+1)+1,sys.iforce1(i)) = dRdlambda( (sys.hforce1(i)-1)*(2*H+1)+1,sys.iforce1(i)) - sys.vforce1(i);
end
for i=1:numel(sys.iforce2)
    dRdlambda( (sys.hforce2(i)-1)*(2*H+1)+1,sys.iforce2(i)) = dRdlambda( (sys.hforce2(i)-1)*(2*H+1)+1,sys.iforce2(i)) - 2*lambda*sys.vforce2(i);
end
for i=1:numel(sys.ic1)
    dRdlambda(1,sys.ic1(i))= dRdlambda(1,sys.ic1(i)) + sys.vc1(i) ;
end
for i=1:numel(sys.ic2)
    dRdlambda(1,sys.ic2(i))= dRdlambda(1,sys.ic2(i)) + sys.vc2(i)*2*lambda ;
end
for i=1:numel(sys.il1)
    dRdlambda(:,sys.il1(i))= dRdlambda(:,sys.il1(i)) + sys.vl1(i)*Ztot(:,sys.jl1(i)) ;
end
dRdlambda=reshape(dRdlambda,coef_per_var*nz_tot,1);




%   dRdomega1
dRdomega1=zeros(coef_per_var,nz_tot);
for i=1:nld
    dRdomega1(:,sys.id(i))= dRdomega1(:,sys.id(i)) + sys.vd(i)*sys.D(Ztot(:,sys.jd(i)),1,0);
end
for i=1:nldl
    dRdomega1(:,sys.idl(i))= dRdomega1(:,sys.idl(i)) + sys.vdl(i)*sys.D(Ztot(:,sys.jdl(i)),1,0);
end
dRdomega1=reshape(dRdomega1,coef_per_var*nz_tot,1);

%   dRdomega2
dRdomega2=zeros(coef_per_var,nz_tot);
for i=1:nld
    dRdomega2(:,sys.id(i))= dRdomega2(:,sys.id(i)) + sys.vd(i)*sys.D(Ztot(:,sys.jd(i)),0,1) ;
end
for i=1:nldl
    dRdomega2(:,sys.idl(i))= dRdomega2(:,sys.idl(i)) + sys.vdl(i)*sys.D(Ztot(:,sys.jdl(i)),0,1);
end
dRdomega2=reshape(dRdomega2,coef_per_var*nz_tot,1);




%   dRdomega1sq
dRdomega1sq =zeros(coef_per_var,nz_tot);
for i=1:nldd
    dRdomega1sq(:,sys.idd(i))= dRdomega1sq(:,sys.idd(i)) + sys.vdd(i)*sys.DD(Ztot(:,sys.jdd(i)),1,0,0) ;
end
dRdomega1sq=reshape(dRdomega1sq,coef_per_var*nz_tot,1);

%   dRdomega2sq
dRdomega2sq =zeros(coef_per_var,nz_tot);
for i=1:nldd
    dRdomega2sq(:,sys.idd(i))= dRdomega2sq(:,sys.idd(i)) + sys.vdd(i)*sys.DD(Ztot(:,sys.jdd(i)),0,1,0) ;
end
dRdomega2sq=reshape(dRdomega2sq,coef_per_var*nz_tot,1);

%   dRdomega1omega2
dRdomega1omega2 =zeros(coef_per_var,nz_tot);
for i=1:nldd
    dRdomega1omega2(:,sys.idd(i))= dRdomega1omega2(:,sys.idd(i)) + sys.vdd(i)*sys.DD(Ztot(:,sys.jdd(i)),0,0,1) ;
end
dRdomega1omega2=reshape(dRdomega1omega2,coef_per_var*nz_tot,1);

%   dRdlambdaomega1
dRdlambdaomega1 =zeros(coef_per_var,nz_tot);
for i=1:nld1
    dRdlambdaomega1(:,sys.id1(i))= dRdlambdaomega1(:,sys.id1(i)) + sys.vd1(i)*sys.D(Ztot(:,sys.jd1(i)),1,0) ;
end
dRdlambdaomega1=reshape(dRdlambdaomega1,coef_per_var*nz_tot,1);

%   dRdlambdaomega2
dRdlambdaomega2 =zeros(coef_per_var,nz_tot);
for i=1:nld1
    dRdlambdaomega2(:,sys.id1(i))= dRdlambdaomega2(:,sys.id1(i)) + sys.vd1(i)*sys.D(Ztot(:,sys.jd1(i)),0,1) ;
end
dRdlambdaomega2=reshape(dRdlambdaomega2,coef_per_var*nz_tot,1);

% Equation defining omega1sq,omega2sq,omega1omega2,lambdaomega1-2
rowomega1sq = sparse([1 1],[sys.neq-1 sys.ninc-4],[-2*omega1 1],1,sys.ninc);
rowomega2sq = sparse([1 1],[sys.neq sys.ninc-3],[-2*omega2 1],1,sys.ninc);
rowomega1omega2 = sparse([1 1 1],[sys.neq-1 sys.neq sys.ninc-2],[-omega2 -omega1 1],1,sys.ninc);
rowlambdaomega1 = sparse([1 1 1],[sys.neq-1 sys.neq+1 sys.ninc-1],[-lambda -omega1 1],1,sys.ninc);
rowlambdaomega2 = sparse([1 1 1],[sys.neq sys.neq+1 sys.ninc],[-lambda -omega2 1],1,sys.ninc);


% phase equation
switch sys.subtype
    case 'autonomous'
        if sys.zi_phase1 <= sys.nz; indice=(sys.zi_phase1-1)*coef_per_var; else, indice=(sys.zi_phase1-1)*coef_per_var+3; end
        rowphase1=sparse(1,indice + (coef_per_var-1)/2+1+(1+2*H),1,1,sys.ninc) ; % sin(omega1) of variable "zi_phase1"
    case 'forced'
        if isfloat(sys.angfreq) % omega1 fixed
            rowphase1=sparse(1,sys.neq-1,1,1,sys.ninc);
        else % omega1 is the continuation parameter
            rowphase1=sparse([1 1],[sys.neq-1 sys.neq+1],[1 -1],1,sys.ninc);
        end
end
if sys.zi_phase2 <= sys.nz; indice=(sys.zi_phase2-1)*coef_per_var; else, indice=(sys.zi_phase2-1)*coef_per_var+3; end
 rowphase2=sparse(1,indice + (coef_per_var-1)/2+2,1,1,sys.ninc) ; % sin(omega2) of variable "zi_phase2"

% Block-decomposition of the Jacobian matrix for the condensation of the
% auxiliar variables

Jacobian.dRdU = [dRdZ(1:sys.neq-2,1:sys.neq-2),sparse(dRdomega1(1:sys.neq-2)),sparse(dRdomega2(1:sys.neq-2)),sparse(dRdlambda(1:sys.neq-2)) ; ...
    rowphase1(1:sys.neq+1) ; rowphase2(1:sys.neq+1)];

Jacobian.dRdUaux = [dRdZ(1:sys.neq-2,sys.neq-1:end), sparse(dRdomega1sq(1:sys.neq-2)) , sparse(dRdomega2sq(1:sys.neq-2)) , ... 
    sparse(dRdomega1omega2(1:sys.neq-2)), sparse(dRdlambdaomega1(1:sys.neq-2)) , sparse(dRdlambdaomega2(1:sys.neq-2)) ; ...
    rowphase1(sys.neq+2:end) ; rowphase2(sys.neq+2:end)];

Jacobian.dRauxdU = [dRdZ(sys.neq-1:end,1:sys.neq-2), sparse(dRdomega1(sys.neq-1:end)) , sparse(dRdomega2(sys.neq-1:end)) , sparse(dRdlambda(sys.neq-1:end)) ; ...
    rowomega1sq(1:sys.neq+1); rowomega2sq(1:sys.neq+1); rowomega1omega2(1:sys.neq+1); rowlambdaomega1(1:sys.neq+1); rowlambdaomega2(1:sys.neq+1)];

Jacobian.dRauxdUaux = [dRdZ(sys.neq-1:end,sys.neq-1:end), sparse(dRdomega1sq(sys.neq-1:end)), sparse(dRdomega2sq(sys.neq-1:end)) , ... 
    sparse(dRdomega1omega2(sys.neq-1:end)), sparse(dRdlambdaomega1(sys.neq-1:end)), sparse(dRdlambdaomega2(sys.neq-1:end)) ; ...
    rowomega1sq(sys.neq+2:end) ;rowomega2sq(sys.neq+2:end) ;rowomega1omega2(sys.neq+2:end) ; rowlambdaomega1(sys.neq+2:end); rowlambdaomega2(sys.neq+2:end)];

Jacobian.dRtotdUtot = [[Jacobian.dRdU ; Jacobian.dRauxdU], [Jacobian.dRdUaux ; Jacobian.dRauxdUaux]];

end


