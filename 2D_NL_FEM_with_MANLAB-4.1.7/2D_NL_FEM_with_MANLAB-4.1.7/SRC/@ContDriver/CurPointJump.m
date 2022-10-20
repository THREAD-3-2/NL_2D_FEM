function [ obj ] = CurPointJump( obj, sys, params )
%Jump from one point to another

disp('Select the jump point');
figure(2); [xjump,yjump] = ginput(1);

U0now=obj.CurPoint.U0now;
Utnow=obj.CurPoint.Utnow;

% if iscell(params.dispvars)
%     switch params.dispvars{params.activecurve,1}
%         case 'lambda'
%             params.dispvars{params.activecurve,1} = sys.neq+1;
%         case 'omega'
%             params.dispvars{params.activecurve,1} = sys.neq;
%     end
%     Znow = sys.get_Ztot(U0now);
%     Ztnow = sys.get_Ztot(Utnow);
%     idxX = params.dispvars{params.activecurve,1}; idxY = params.dispvars{params.activecurve,2};
%     xdep = U0now(idxX); ydep = sqrt(Znow(:,idxY)'*Znow(:,idxY));
%     txdep = Utnow(idxX); tydep = Znow(:,idxY)'*Ztnow(:,idxY)/ydep;
%     % Normalisation and orientation of the printed tangent
%     Utnow = Utnow/sqrt((txdep^2+tydep^2));
%     Ztnow = sys.get_Ztot(Utnow);
%     txnow = Utnow(idxX); tynow = Znow(:,idxY)'*Ztnow(:,idxY)/ydep;
%     scal = (xjump-xdep)*txnow + (yjump-ydep)*tynow;
% else
idx = params.dispvars(params.activecurve,1);
idy = params.dispvars(params.activecurve,2);
xdep = U0now(idx); ydep = U0now(idy);
% Normalisation and orientation of the printed tangent
Utnow = Utnow/sqrt((Utnow(idx)^2+Utnow(idy)^2));
scal = (xjump-xdep)*Utnow(idx) + (yjump-ydep)*Utnow(idy);
%end

if scal < 0, Utnow = -Utnow; end
% Effective jump
jump_length = sqrt((xjump-xdep)^2+(yjump-ydep)^2);
U0now = U0now + jump_length*Utnow;

disp('Correction to go back on the branch of solution...');
U0now = NRcorrection(sys, U0now);
Utnow = tangentvector(sys,Jacobian(sys,U0now));

obj.CurPoint.U0now=U0now;
obj.CurPoint.Utnow=Utnow;
obj.CurPoint.Status='regular';
