function Jacobian = Jacobian_HBM(sys,Utot)
% Compute the Jacobian of R with respect to U

H    = sys.H;
nz_tot   = sys.nz_tot;

DHp1=2*H+1 ;
DHp12=DHp1*DHp1;


% Extract informations from Utot
[Ztot,omega,lambda,omega2,lambdaomega] = sys.get_Ztot(Utot);

Ztot_t0 = sum(Ztot(1:H+1,:),1);

% definition of the small (constant)matrix : I and d
%
% Identity matrix (2*H+1) sparse
[iI,jI,vI]=find(speye(DHp1)) ;

% d sparse matrix
idmat=[ 1 ,   2:H+1  ,H+2:DHp1 ]';
jdmat=[ 1 , H+2:DHp1 ,  2:H+1  ]';
vdmat=[ 0 ,     1:H  , -(1:H)  ]';

% dd sparse matrix
iddmat=[ 1 ,   2:H+1  ,H+2:DHp1 ]';
jddmat=iddmat;
vddmat=-[ 0 ,     (1:H).^2  , (1:H).^2  ]';

% ib and jb for b sparce matrix
[ib,jb,~]=find( sparse( ones(DHp1,DHp1) ));

%   dRdZ =  L0 + lambda L1  + 2*Q(Z,.)
%   construct the list for the large matrix L0, L1, dl, 2*Q(Z,.) and dq(Z,.)
nl0=numel(sys.il0); iL0=zeros(DHp1,nl0); jL0=zeros(DHp1,nl0); vL0=zeros(DHp1,nl0);
for i=1:nl0
    iL0(:,i)=iI + (sys.il0(i)-1)*DHp1 ;
    jL0(:,i)=jI + (sys.jl0(i)-1)*DHp1 ;
    vL0(:,i)=vI*sys.vl0(i)  ;
end

nl1=numel(sys.il1); iL1=zeros(DHp1,nl1); jL1=zeros(DHp1,nl1); vL1=zeros(DHp1,nl1);
for i=1:nl1
    iL1(:,i)=iI + (sys.il1(i)-1)*DHp1 ;
    jL1(:,i)=jI + (sys.jl1(i)-1)*DHp1 ;
    vL1(:,i)=vI*(sys.vl1(i)*lambda) ;
end

nq=numel(sys.iq); iQL=zeros(DHp12,nq); jQL=zeros(DHp12,nq); vQL=zeros(DHp12,nq);
for i=1:nq
    % Avoid to compute the same B several times / it uses that obj.jq is non-decreasing.
    if ~(i>1 && sys.jq(i-1) == sys.jq(i))
        vb= reshape(sys.B(Ztot(:,sys.jq(i))) ,DHp12,1) ;
    end
    iQL(:,i)=ib + (sys.iq(i)-1)*DHp1 ;
    jQL(:,i)=jb + (sys.kq(i)-1)*DHp1 ;
    vQL(:,i)=2*vb*sys.vq(i) ;
end

nld=numel(sys.id); iD=zeros(DHp1,nld); jD=zeros(DHp1,nld); vD=zeros(DHp1,nld);
for i=1:nld
    iD(:,i)=idmat + (sys.id(i)-1)*DHp1 ;
    jD(:,i)=jdmat + (sys.jd(i)-1)*DHp1 ;
    vD(:,i)=vdmat*(sys.vd(i)*omega) ;
end

nldd=numel(sys.idd); idd=zeros(DHp1,nldd); jdd=zeros(DHp1,nldd); vdd=zeros(DHp1,nldd);
for i=1:nldd
    idd(:,i)=iddmat + (sys.idd(i)-1)*DHp1 ;
    jdd(:,i)=jddmat + (sys.jdd(i)-1)*DHp1 ;
    vdd(:,i)=vddmat*(sys.vdd(i)*omega2) ;
end

ndl=numel(sys.idl); idl=zeros(DHp1,ndl); jdl=zeros(DHp1,ndl); vdl=zeros(DHp1,ndl);
for i=1:ndl
    idl(:,i)=idmat + (sys.idl(i)-1)*DHp1 ;
    jdl(:,i)=jdmat + (sys.jdl(i)-1)*DHp1 ;
    vdl(:,i)=vdmat*sys.vdl(i) ;
end

idq1mat = repmat(idmat,DHp1,1);
jdq1mat = reshape(repmat(jdmat',DHp1,1),DHp12,1);
vdq1mat = reshape(repmat(vdmat',DHp1,1),DHp12,1);
ndq=numel(sys.idq); idq1=zeros(DHp12,ndq); jdq1=zeros(DHp12,ndq); vdq1=zeros(DHp12,ndq);
for i=1:ndq
    vb= reshape(sys.B(Ztot(:,sys.jdq(i))),DHp12,1) ;
    idq1(:,i)=idq1mat + (sys.idq(i)-1)*DHp1 ;
    jdq1(:,i)=jdq1mat + (sys.kdq(i)-1)*DHp1 ;
    vdq1(:,i)=(vdq1mat.*vb)*sys.vdq(i) ;
end

idq2=zeros(DHp12,ndq); jdq2=zeros(DHp12,ndq); vdq2=zeros(DHp12,ndq);
for i=1:ndq
    vb= reshape(sys.B(sys.D(Ztot(:,sys.kdq(i)))),DHp12,1) ;
    idq2(:,i)=ib + (sys.idq(i)-1)*DHp1 ;
    jdq2(:,i)=jb + (sys.jdq(i)-1)*DHp1 ;
    vdq2(:,i)=vb*sys.vdq(i) ;
end

nld1=numel(sys.id1); id1=zeros(DHp1,nld1); jd1=zeros(DHp1,nld1); vd1=zeros(DHp1,nld1);
for i=1:nld1
    id1(:,i)=idmat + (sys.id1(i)-1)*DHp1 ;
    jd1(:,i)=jdmat + (sys.jd1(i)-1)*DHp1 ;
    vd1(:,i)=vdmat*(sys.vd1(i)*lambdaomega) ;
end

idRdZ=[ reshape(iL0,nl0*DHp1,1) ; reshape(iL1,nl1*DHp1,1) ; ...
    reshape(iQL,nq*DHp12,1) ; reshape(iD,nld*DHp1,1); ...
    reshape(id1,nld1*DHp1,1) ; reshape(idd,nldd*DHp1,1); ...
    reshape(idq2,ndq*DHp12,1) ; reshape(idq1,ndq*DHp12,1) ; reshape(idl,ndl*(DHp1),1)];

jdRdZ=[ reshape(jL0,nl0*DHp1,1) ; reshape(jL1,nl1*DHp1,1) ; ...
    reshape(jQL,nq*DHp12,1) ; reshape(jD,nld*DHp1,1); ...
    reshape(jd1,nld1*DHp1,1) ; reshape(jdd,nldd*DHp1,1); ...
    reshape(jdq2,ndq*DHp12,1) ; reshape(jdq1,ndq*DHp12,1) ; reshape(jdl,ndl*(DHp1),1)];

vdRdZ=[ reshape(vL0,nl0*DHp1,1) ; reshape(vL1,nl1*DHp1,1) ; ...
    reshape(vQL,nq*DHp12,1) ; reshape(vD,nld*DHp1,1); ...
    reshape(vd1,nld1*DHp1,1) ; reshape(vdd,nldd*DHp1,1); ...
    reshape(vdq2,ndq*DHp12,1) ; reshape(vdq1,ndq*DHp12,1) ; reshape(vdl,ndl*(DHp1),1)];

dRdZ=sparse(idRdZ,jdRdZ,vdRdZ,DHp1*nz_tot,DHp1*nz_tot) ;


%  Generalization with mass matrix

iM=zeros(DHp1,nld); jM=zeros(DHp1,nld); vM=zeros(DHp1,nld);
for i=1:nld
    iM(:,i)=iI + (sys.id(i)-1)*DHp1 ;
    jM(:,i)=jI + (sys.jd(i)-1)*DHp1 ;
    vM(:,i)=vI*(sys.vd(i)) ;
end

Mfull=sparse(reshape(iM,nld*DHp1,1),reshape(jM,nld*DHp1,1),reshape(-vM,nld*DHp1,1),DHp1*nz_tot,DHp1*nz_tot);


%   Modification for dl/dq : initial condition replace the zero-th harmonic :
dRdZ(DHp1*(sys.idl-1)+1,:) = 0;
for i=1:numel(sys.idl)
    dRdZ(DHp1*(sys.idl(i)-1)+1,DHp1*(sys.jdl(i)-1)+(1:H+1)) = sys.vdl(i);
end
for i=1:numel(sys.idq)
    dRdZ(DHp1*(sys.idq(i)-1)+1,DHp1*(sys.kdq(i)-1)+(1:H+1)) = sys.vdq(i)*Ztot_t0(sys.jdq(i))  +  dRdZ(DHp1*(sys.idq(i)-1)+1,DHp1*(sys.kdq(i)-1)+(1:H+1));
end

%   dRdomega  = D(Z)
dRdomega =zeros(DHp1,nz_tot);
for i=1:nld
    dRdomega(:,sys.id(i))= dRdomega(:,sys.id(i)) + sys.vd(i)*sys.D(Ztot(:,sys.jd(i))) ;
end
dRdomega=reshape(dRdomega,DHp1*nz_tot,1);

%   dRdomega2  = dd(Z)
dRdomega2 =zeros(DHp1,nz_tot);
for i=1:nldd
    dRdomega2(:,sys.idd(i))= dRdomega2(:,sys.idd(i)) + sys.vdd(i)*sys.D(sys.D(Ztot(:,sys.jdd(i)))) ;
end
dRdomega2=reshape(dRdomega2,DHp1*nz_tot,1);
% equation defining omega2
rowomega2 = sparse([1 1],[sys.neq sys.ninc-1],[-2*omega 1],1,sys.ninc);

%   dRd1ambda  = C1 + L1(Z) + 2*lambda*C2
dRdlambda =zeros(DHp1,nz_tot);
for i=1:numel(sys.iforce1)
    dRdlambda(sys.hforce1(i),sys.iforce1(i)) = dRdlambda(sys.hforce1(i),sys.iforce1(i)) - sys.vforce1(i);
end
for i=1:numel(sys.ic1)
    dRdlambda(1,sys.ic1(i))= dRdlambda(1,sys.ic1(i)) + sys.vc1(i) ;
end
for i=1:numel(sys.iforce2)
    dRdlambda(sys.hforce2(i),sys.iforce2(i)) = dRdlambda(sys.hforce2(i),sys.iforce2(i)) - sys.vforce2(i)*2*lambda;
end
for i=1:numel(sys.ic2)
    dRdlambda(1,sys.ic2(i))= dRdlambda(1,sys.ic2(i)) + sys.vc2(i)*2*lambda ;
end
for i=1:numel(sys.il1)
    dRdlambda(:,sys.il1(i))= dRdlambda(:,sys.il1(i)) + sys.vl1(i)*Ztot(:,sys.jl1(i)) ;
end
dRdlambda=reshape(dRdlambda,DHp1*nz_tot,1);

%   dRd1ambdaomega  = d1(Z)
dRdlambdaomega =zeros(DHp1,nz_tot);
for i=1:numel(sys.id1)
    dRdlambdaomega(:,sys.id1(i))= dRdlambdaomega(:,sys.id1(i)) + sys.vd1(i)*sys.D(Ztot(:,sys.jd1(i))) ;
end
dRdlambdaomega=reshape(dRdlambdaomega,DHp1*nz_tot,1);
% Equation defining lambdaomega
rowlambdaomega = sparse([1 1 1],[sys.neq sys.neq+1 sys.ninc],[-lambda -omega 1],1,sys.ninc);

% phase equation zi'(0)=0 -> zis1 + 2*zis2 + 3*zis3 + ...  =0
switch sys.subtype
    case 'autonomous'
        if sys.zi_phase1 <= sys.nz; indice=(sys.zi_phase1-1)*DHp1; else, indice=(sys.zi_phase1-1)*DHp1+2; end
        %%% derivative of zi_phase1 variable null at t=0.
        rowphase=sparse( ones(1,H) ,indice+H+2:indice+DHp1,(1:H),1,sys.ninc) ;
        %%% first sine of zi_phase1 variable null.
        %rowphase=sparse(1,indice+H+2,1,1,obj.ninc) ;
    case 'forced'
        if isfloat(sys.angfreq) % omega fixed
            rowphase=sparse(1,sys.neq,1,1,sys.ninc);
        else % omega is the continuation parameter
            rowphase=sparse( ones(1,2) ,sys.neq:sys.neq+1,[1 -1],1,sys.ninc);
        end
end

% Block-decomposition of the Jacobian matrix for the condensation of the
% auxiliar variables

Jacobian.dRdU = [dRdZ(1:sys.neq-1,1:sys.neq-1),sparse(dRdomega(1:sys.neq-1)),sparse(dRdlambda(1:sys.neq-1)) ; rowphase(1:sys.neq+1)];

Jacobian.dRdUaux = [dRdZ(1:sys.neq-1,sys.neq:end), sparse(dRdomega2(1:sys.neq-1)) , sparse(dRdlambdaomega(1:sys.neq-1)) ; rowphase(sys.neq+2:end)];
Jacobian.dRauxdU = [dRdZ(sys.neq:end,1:sys.neq-1), sparse(dRdomega(sys.neq:end)) , sparse(dRdlambda(sys.neq:end)) ; rowomega2(1:sys.neq+1) ; rowlambdaomega(1:sys.neq+1)];
Jacobian.dRauxdUaux = [dRdZ(sys.neq:end,sys.neq:end), sparse(dRdomega2(sys.neq:end)) , sparse(dRdlambdaomega(sys.neq:end)) ; rowomega2(sys.neq+2:end) ; rowlambdaomega(sys.neq+2:end)];

Jacobian.dRtotdUtot = [[Jacobian.dRdU ; Jacobian.dRauxdU], [Jacobian.dRdUaux ; Jacobian.dRauxdUaux]];

% Generalization with mass matrix

Jacobian.Mfull = Mfull;
Jacobian.M = [Mfull(1:sys.neq-1,1:sys.neq-1), sparse([],[],[],sys.neq-1,2); sparse([],[],[],1,sys.neq+1)] ;
Jacobian.Maux = [Mfull(sys.neq:end,1:sys.neq-1),sparse([],[],[],sys.neq_aux-2,2);sparse([],[],[],2,sys.neq+1)];


end


