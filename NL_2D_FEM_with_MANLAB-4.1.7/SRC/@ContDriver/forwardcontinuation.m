function [ obj ] = forwardcontinuation( obj, sys, params, nbforward )
%   performs nbforwrd continuation steps, creating nbforward CheckPoint object

U0now=obj.CurPoint.U0now ;     % get the ''current point'' of the diagram object
Utnow=obj.CurPoint.Utnow  ;
BifStatus=obj.CurPoint.BifStatus;

for i=1:nbforward
    
    if get(sys,'NRmethod') ~= 0
        U0now = NRcorrection(sys,U0now);
    end
    
    switch BifStatus
        case 'simplebif'
            [U0,Amax,DataBif,DataStab]    = ANMserieBif(sys,U0now,Utnow,obj.CurPoint.Ut2);
        otherwise
            [U0,Amax,DataBif,DataStab]    = ANMseries(sys,U0now,Utnow)  ;    %  U0 (Taylor type)
    end
    if get(sys,'BifDetection') == 0
        DataBif.status = 'nothing';
    end
    ChckPoint      = CheckPoint(U0,Utnow,Amax,params,DataBif,DataStab);
    
    R0_tot = sys.R(sys,U0now) + sys.PerturbationSize*sys.pertvect;
    Rend_tot = sys.R(sys,ChckPoint.Uend) + sys.PerturbationSize*sys.pertvect;
    
    disp(['Section ' num2str(i) ': ||R(U(a=0))||=' num2str(norm(R0_tot))  ...
        '||R(U(a=Amax))||=' num2str(norm(Rend_tot)) ';  Amax=' num2str(Amax)]);
    
    obj=CheckPointAdd(obj,ChckPoint);
    
    if params.alldisplayoff == 0
        
        disp(ChckPoint);
        if params.globaldisplay == 1
            sys.global_display(sys,ChckPoint);
        end
        drawnow;
        
    end
    
    U0now=ChckPoint.Uend;
    Utnow=ChckPoint.Utend;
    BifStatus='regular';
end

if params.pointdisplay == 1
    sys.point_display(sys, U0now); drawnow;
end

obj.CurPoint.U0now=U0now;   % set the last end point  in the ''current point''
obj.CurPoint.Utnow=Utnow;
obj.CurPoint.BifStatus='regular';

