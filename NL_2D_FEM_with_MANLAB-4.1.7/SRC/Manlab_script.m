function Diagram = Manlab_script(varargin)

pt_disp = 0;
gb_disp= 0;
nbforward_step=10;

propertyArgin = varargin;
while length(propertyArgin) >=2
    prop = propertyArgin{1};
    val=propertyArgin{2};
    propertyArgin=propertyArgin(3:end);
    switch prop
        case 'sys'
            sys=val;
        case 'type'
            type=val;
        case 'U0value'
            U0value=val;
            [m,n]=size(U0value);
            if (n>1)
                errordlg('Dimension mismatch: U0value should be a column vector.','on');
            end
        case 'nb_step'
            nbforward_step=val;
        case 'displayvariables'
            dispvars=val;
        case 'name'
            sys=set(sys,'name',val);
        case 'order'
            sys=set(sys,'order',val);
        case 'ANMthreshold'
            sys=set(sys,'ANMthreshold',val);
        case 'NRthreshold'
            sys=set(sys,'NRthreshold',val);
        case 'NRmethod'
            sys=set(sys,'NRmethod',val);
        case 'BifDetection'
            sys=set(sys,'BifDetection',val);
        case 'StabilityCheck'
            sys=set(sys,'StabilityCheck',val);
        case 'StabTol'
            sys=set(sys,'StabTol',val);
        case 'NRitemax'
            sys=set(sys,'NRitemax',val);
        case 'Amax_max'
            sys=set(sys,'Amax_max',val);
        case 'PointDisplay'
            pt_disp = val;
        case 'GlobalDisplay'
            gb_disp = val;
    end
end

% CURRENT POINT
U0now = NRcorrection(sys,U0value);       % Perform correction on U0value if needed
Jacobian_structure = sys.Jacobian(U0now);
[Utnow,~,Jacobian_structure] = sys.tangentvector(Jacobian_structure);   % now compute the tangent at U0now
if sys.StabilityCheck == 1
    stab_flag = sys.Stability(U0value,Jacobian_structure);
    if stab_flag == 1; CurPoint.StabStatus = 'stable';
    else; CurPoint.StabStatus = 'unstable'; end
end
CurPoint.U0now = U0now    ;        %  set U0now in the Current point
CurPoint.Utnow = Utnow    ;        %  set Utnow
CurPoint.Ut2   = [ ]        ;              %  For storing a second tangent when the CurPoint is a simple bif
CurPoint.BifStatus= 'regular'  ;              % 'regular' or 'simplebif'

% DIAGRAM
diagram_struct = ContDriver(CurPoint);

% PLOT PARAMETRES
params.axisauto   = 1;
params.uaxis      = [0 1 0 1];
params.drawtype   = {'-',':'};
params.dispcolors = ['b','r','g','c','m','y','k'];
params.nbpts      = 10;
params.dispvars   = dispvars;
params.markers    = '.';
params.globaldisplay   = gb_disp;
params.pointdisplay    = pt_disp;

params.alldisplayoff = 1;

% Display diagram
fig2=figure(2);
disp(diagram_struct);

hcurpoint=CurPointDisp(diagram_struct,params,sys);

scrsz = get(0,'ScreenSize');
left=scrsz(3)/7;bottom=scrsz(4)/3; width=scrsz(3)*2/5; height=scrsz(4)/2;
set( fig2 , 'position' , [left bottom width height]);

fprintf('Action. Forward computation of %i sections \n',nbforward_step);

CurPointErase(diagram_struct,params,hcurpoint);

diagram_struct=forwardcontinuation(diagram_struct,sys,params,nbforward_step);

if pt_disp == 1
    CurPointDisp(diagram_struct,params,sys);
end

Diagram = diagram_struct.ChckPoint;

end