function  Manlab(varargin)
% MANLAB MATLAB code for Manlab.fig
%      MANLAB, by itself, creates a new MANLAB or raises the existing
%      singleton*.
%
%      H = MANLAB returns the handle to a new MANLAB or the handle to
%      the existing singleton*.
%
%      MANLAB('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MANLAB.M with the given input arguments.
%
%      MANLAB('Property','Value',...) creates a new MANLAB or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Manlab_openingfcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Manlab_openingfcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Manlab

% Last Modified by GUIDE v2.5 24-Mar-2019 10:24:23

% Begin initialization code - DO NOT EDIT

fprintf('\n');
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @Manlab_OpeningFcn, ...
    'gui_OutputFcn',  @Manlab_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

% --- Executes just before Manlab is made visible.
function Manlab_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Manlab (see VARARGIN)

% Default display.
pt_disp = 1;
gb_disp = 1;
nbpts = 15;

% Default initial correction.
firstcorr = 1;

propertyArgin = varargin;
while length(propertyArgin) >=2
    prop = propertyArgin{1};
    val=propertyArgin{2};
    propertyArgin=propertyArgin(3:end);
    switch prop
        case 'sys'
            sys=val;
        case 'U0value'
            U0value=val;
            [m,n]=size(U0value);
            if (n>1)
                errordlg('Dimension mismatch: U0value should be a column vector.','on');
            end
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
        case 'NRstart'
            firstcorr = val;
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
        case 'nbpts'
            nbpts = val;
    end
end

% SYSTEM
handles.sys=sys;

% CURRENT POINT
if firstcorr == 1
    handles.U0now = NRcorrection(handles.sys,U0value);       % Perform correction on U0value if needed
else
    handles.U0now = U0value;
end

Jacobian_structure = handles.sys.Jacobian(handles.U0now);
[handles.Utnow,~,Jacobian_structure] = handles.sys.tangentvector(Jacobian_structure);   % now compute the tangent at U0now

% J = handles.sys.Jacobian(handles.U0now).dRtotdUtot;
% en = zeros(1, size(J,2));
% en((2*sys.H + 1)*sys.nz+1) = 1;
% Jac = [ J;
%         en];
% J(:,(2*sys.H + 1)*sys.nz+1)=[];
% Jacobian = sys.Jacobian(handles.U0now);
% [Utf,Vtf,Jacobian] = handles.sys.tangentvector(Jacobian);
if isfield(sys.parameters,'model')    
[Z0,omega,lambda,omega2,lambdaomega] = get_Ztot(sys,handles.U0now);
[Zt,omegat,lambdat,omega2t,lambdaomegat] = get_Ztot(sys,handles.Utnow);
q0 = zeros(3*handles.sys.parameters.model.mesh.number_nodes,3);
qt = zeros(3*handles.sys.parameters.model.mesh.number_nodes,3);
active_dof = handles.sys.parameters.model.boundary.active_dof;
n_act = length(active_dof);
q0(active_dof,1) = Z0(1,1:n_act)'; % constant 
q0(active_dof,2) = Z0(2,1:n_act)'; % harmonic 1 cosine
q0(active_dof,3) = Z0(2+sys.H,1:n_act)';% harmonic 1 sine
qt(active_dof,1) = Zt(1,1:n_act)'; % constant 
qt(active_dof,2) = Zt(2,1:n_act)'; % harmonic 1 cosine
qt(active_dof,3) = Zt(2+sys.H,1:n_act)';% harmonic 1 sine
fig0 = figure(88);
subplot(2,1,1); hold on
handles.sys.parameters.model.plot_deformed_mesh(q0(:,1),fig0,'-k')
handles.sys.parameters.model.plot_deformed_mesh(q0(:,2),fig0,'-r')
handles.sys.parameters.model.plot_deformed_mesh(q0(:,3),fig0,'-b')
legend('H0', 'H1', 'H2')
ylabel('shape at initial point')
subplot(2,1,2); hold on
handles.sys.parameters.model.plot_deformed_mesh(qt(:,1),fig0,'-k')
handles.sys.parameters.model.plot_deformed_mesh(qt(:,2),fig0,'-r')
handles.sys.parameters.model.plot_deformed_mesh(qt(:,3),fig0,'-b')
legend('H0', 'H1', 'H2')
ylabel('tangent vector at initial point')
end

if sys.StabilityCheck == 1
    stab_flag = handles.sys.Stability(U0value,Jacobian_structure);
    if stab_flag == 1; CurPoint.StabStatus = 'stable';
    else; CurPoint.StabStatus = 'unstable'; end
end
CurPoint.U0now = handles.U0now    ;        %  set U0now in the Current point
CurPoint.Utnow = handles.Utnow    ;        %  set Utnow
CurPoint.Ut2   = [ ]        ;              %  For storing a second tangent when the CurPoint is a simple bif
CurPoint.BifStatus= 'regular'  ;              % 'regular' or 'simplebif'
handles.CurPoint=CurPoint ;
handles.CurPointInit=CurPoint ;

% DIAGRAM
handles.diagram = ContDriver(handles.CurPoint);


% PLOT PARAMETRES
params.axisauto   = 1;
params.uaxis      = [0 1 0 1];
params.drawtype   = {'-',':'}; % {'--','-.'} for perturbed branches, see EdPerturbationSize below.
params.dispcolors = ['b','r','g','c','m','y','k','b'];
params.nbpts      = nbpts;
params.dispvars   = dispvars;

params.alldisplayoff = 0;

handles.params=params;



% Display diagram
fig2=figure(2);
disp(handles.diagram);
handles.hcurpoint=CurPointDisp(handles.diagram,handles.params,handles.sys);

scrsz = get(0,'ScreenSize');
left=scrsz(3)/7;bottom=scrsz(4)/3; width=scrsz(3)*2/5; height=scrsz(4)/2;
figure('MenuBar','none','Name','Manlab','NumberTitle','off',...
    'visible','off');
set( fig2 , 'position' , [left bottom width height]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Continuation management through the interface
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data for the GUI at start
handles.nbforward    = 5;
handles.NRmethod     = get(handles.sys,'NRmethod');
handles.NRthreshold  = get(handles.sys,'NRthreshold');
handles.ANMthreshold = get(handles.sys,'ANMthreshold');
handles.Amax_max = get(handles.sys,'Amax_max');
handles.PerturbationSize = get(handles.sys,'PerturbationSize');
handles.BifDetection = get(handles.sys,'BifDetection');
handles.StabilityCheck = get(handles.sys,'StabilityCheck');
handles.CancelLast = 1;
handles.params.activecurve  = 1;
handles.params.markers      = 1;
handles.params.lockzoom     = 0;
handles.params.pointdisplay  = pt_disp;
handles.params.globaldisplay  = gb_disp;

%Setting the values on the UI editbox.

set(handles.EdANMthreshold,'string',num2str(handles.sys.ANMthreshold));
set(handles.EdAmax_max,'string',num2str(handles.sys.Amax_max));
set(handles.EdCorrection,'string',num2str(handles.sys.NRthreshold));
set(handles.EdCorrType,'string',num2str(handles.sys.NRmethod));
set(handles.EdPerturbationSize,'string',num2str(handles.sys.PerturbationSize));
set(handles.EdCancelLast,'string',num2str(handles.CancelLast));
set(handles.EdSection,'string',handles.nbforward);
set(handles.EdJump,'string', num2str(handles.params.activecurve));

%Setting Bifurcation check box.
if handles.sys.BifDetection==0
    set(handles.Bifurcation, 'Value',0)
else
    set(handles.Bifurcation, 'value',1)
end

%Setting Stability check box.
if handles.sys.StabilityCheck==0
    set(handles.Stability, 'Value',0)
else
    set(handles.Stability, 'value',1)
end

%Setting point display check box.
if handles.params.pointdisplay==0
    set(handles.CkbPointdisplay, 'Value',0)
else
    set(handles.CkbPointdisplay, 'value',1)
end

%Setting global display check box.
if handles.params.globaldisplay==0
    set(handles.CkbGlobaldisplay, 'Value',0)
else
    set(handles.CkbGlobaldisplay, 'value',1)
end

% Setting curve markers check box.
if handles.params.markers==0
    set(handles.CkbMarkers, 'Value',0)
else
    set(handles.CkbMarkers, 'value',1)
end


% Update handles structure
% Choose default command line output for Manlab
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Manlab wait for user response (see UIRESUME)
% uiwait(handles.Manlab);


% --- Outputs from this function are returned to the command line.
function varargout = Manlab_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


function BtCancelSection_Callback(hObject, eventdata, handles)
% hObject    handle to BtCancelSection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.diagram = CheckPointCancel(handles.diagram,'one');
fprintf('Action. Cancel section ');
disp(handles.diagram);
handles.hcurpoint=CurPointDisp(handles.diagram,handles.params,handles.sys);
guidata(hObject, handles);

function BtCancelLast_Callback(hObject, eventdata, handles)
% hObject    handle to BtCancelSection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.diagram = CheckPointCancel(handles.diagram,'last',handles.CancelLast);
fprintf(['Action. Cancel last ' num2str(handles.CancelLast) ' section(s).']);
disp(handles.diagram);
handles.hcurpoint=CurPointDisp(handles.diagram,handles.params,handles.sys);
guidata(hObject, handles);

% --- Executes on button press in BtCancelAll.
function BtCancelAll_Callback(hObject, eventdata, handles)
% hObject    handle to BtCancelAll (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.diagram = CheckPointCancel(handles.diagram,'all');
fprintf('Action. Cancel the diagram');
disp(handles.diagram);
handles.hcurpoint=CurPointDisp(handles.diagram,handles.params,handles.sys);
guidata(hObject, handles);

% --- Executes on button press in BtLoad.
function BtLoad_Callback(hObject, eventdata, handles)
% hObject    handle to BtLoad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[file,path] = uigetfile;
if (~isempty(file) && ~isempty(path))
    loadeddatas = load('-mat',[path file],'diagram') ;
    handles.diagram = loadeddatas.diagram;
    disp(handles.diagram);
    fprintf('Action. Load  diagram');
end
guidata(hObject, handles);

% --- Executes on button press in BtSave.
function BtSave_Callback(hObject, eventdata, handles)
% hObject    handle to BtSave (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[file,path] = uiputfile;
diagram = handles.diagram;
save([path file], 'diagram');
fprintf('Action. Save  diagram');
guidata(hObject, handles);

% --- Executes on button press in BtVariables.
function BtVariables_Callback(hObject, eventdata, handles)
% hObject    handle to BtVariables (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.hcurpoint=CurPointErase( handles.diagram, handles.params, handles.hcurpoint);
[handles.diagram,handles.params]=displayvariables(handles.diagram,handles.params);
handles.hcurpoint= CurPointDisp( handles.diagram, handles.params);
fprintf('Action. change display variables and redraw  ');

guidata(hObject, handles);

% --- Executes on button press in CkbMarkers.
function CkbMarkers_Callback(hObject, eventdata, handles)
% hObject    handle to CkbMarkers (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.params.markers=mod(handles.params.markers+1,2);
if handles.params.markers==1
    fprintf('Action. Markers On')
else
    fprintf('Action. Markers Off')
end
guidata(hObject, handles);
% Hint: get(hObject,'Value') returns toggle state of CkbMarkers


% --- Executes on button press in CkbZoom.
function CkbZoom_Callback(hObject, eventdata, handles)
% hObject    handle to CkbZoom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of CkbZoom


% --- Executes on button press in BtSelectPoint.
function BtSelectPoint_Callback(hObject, eventdata, handles)
% hObject    handle to BtSelectPoint (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
fprintf('Action: Display of the selected point');
handles.U = visupoint(handles.diagram,handles.params);
handles.sys.point_display(handles.sys, handles.U ); drawnow;
guidata(hObject, handles);
% --- Executes on button press in BtSelectPoint.

function BtInitial_Callback(hObject, eventdata, handles)
% hObject    handle to BtInitial (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.hcurpoint=CurPointErase( handles.diagram, handles.params, handles.hcurpoint);
handles.diagram.CurPoint=handles.CurPointInit;
handles.hcurpoint= CurPointDisp( handles.diagram, handles.params ,handles.sys);
if  handles.params.pointdisplay ==1
    handles.sys.point_display( handles.sys, get(handles.diagram,'U0now'));  drawnow;
end
fprintf('Action. replace current point to initial position');
guidata(hObject, handles);
% --- Executes on button press in BtInitial.

function BtSet_Callback(hObject, eventdata, handles)
% hObject    handle to BtSet (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
fprintf('Action. set the current point \n');
figure(2);
if length(get(handles.diagram,'ChckPoint')) == 0
    disp('Empty diagram');
else
    handles.hcurpoint=CurPointErase(handles.diagram,handles.params,handles.hcurpoint);
    [handles.diagram,handles.params,handles.sys]=CurPointSet(handles.diagram,handles.params,handles.sys);
    
    handles.hcurpoint=CurPointDisp(handles.diagram,handles.params,handles.sys);
    if handles.params.pointdisplay ==1
        handles.sys.point_display(handles.sys, get(handles.diagram,'U0now')); drawnow;
    end
end
guidata(hObject, handles);

function BtImport_Callback(hObject, eventdata, handles)
% hObject    handle to BtImport (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global U
handles.hcurpoint=CurPointErase( handles.diagram, handles.params, handles.hcurpoint);
handles.diagram.CurPoint.U0now = U;
handles.diagram.CurPoint.Utnow = handles.sys.tangentvector(handles.sys.Jacobian(U));
handles.diagram.CurPoint.Status='regular';
handles.hcurpoint= CurPointDisp( handles.diagram, handles.params);
if  handles.params.pointdisplay ==1
    handles.sys.point_display( handles.sys, get(handles.diagram,'U0now'));  drawnow;
end
fprintf('Action. replace current point by point in global variable U');
guidata(hObject, handles);
% --- Executes on button press in BtImport.


% --- Executes on button press in Btjump.
function Btjump_Callback(hObject, eventdata, handles)
% hObject    handle to Btjump (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.diagram=CurPointJump(handles.diagram,handles.sys,handles.params);
handles.hcurpoint=CurPointErase(handles.diagram,handles.params,handles.hcurpoint);
fprintf('Action: Jump to the selected point');
handles.hcurpoint=CurPointDisp(handles.diagram,handles.params,handles.sys);
if handles.params.pointdisplay ==1
    handles.sys.point_display(handles.sys, get(handles.diagram,'U0now')); drawnow;
end
guidata(hObject, handles);


% --- Executes on button press in BtInvTangent.
function BtInvTangent_Callback(hObject, eventdata, handles)
% hObject    handle to BtInvTangent (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.hcurpoint=CurPointErase(handles.diagram,handles.params,handles.hcurpoint);
handles.diagram=CurPointInvTan(handles.diagram);
fprintf('Action. Invert tangent vector')
handles.hcurpoint=CurPointDisp(handles.diagram,handles.params,handles.sys);
guidata(hObject, handles);

function EdSection_Callback(hObject, eventdata, handles)
% hObject    handle to EdSection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of EdSection as text
%        str2double(get(hObject,'String')) returns contents of EdSection as a double
handles.nbforward=str2num(get(handles.EdSection,'String'));
fprintf('Action: Number of forward computations updated to %i ',handles.nbforward);
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function EdSection_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EdSection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in BtForward.
function BtForward_Callback(hObject, eventdata, handles)
% hObject    handle to BtForward (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

fprintf('Action. Forward computation of %i sections \n',handles.nbforward);
handles.hcurpoint=CurPointErase(handles.diagram,handles.params,handles.hcurpoint);
handles.diagram=forwardcontinuation(handles.diagram,handles.sys,handles.params,handles.nbforward);
handles.hcurpoint=CurPointDisp(handles.diagram,handles.params,handles.sys);

guidata(hObject, handles);


function EdANMthreshold_Callback(hObject, eventdata, handles)
% hObject    handle to EdANMthreshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.sys.ANMthreshold=str2double(get(handles.EdANMthreshold,'String'));
fprintf('ANM threshold is updated to %i ',handles.sys.ANMthreshold);
guidata(hObject, handles);
% Hints: get(hObject,'String') returns contents of EdANMthreshold as text
%        str2double(get(hObject,'String')) returns contents of EdANMthreshold as a double


% --- Executes during object creation, after setting all properties.
function EdANMthreshold_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EdANMthreshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function EdAmax_max_Callback(hObject, eventdata, handles)
% hObject    handle to EdANMthreshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.sys.Amax_max=str2double(get(handles.EdAmax_max,'String'));
fprintf('Maximum value of Amax is updated to %i ',handles.sys.Amax_max);
guidata(hObject, handles);
% Hints: get(hObject,'String') returns contents of EdANMthreshold as text
%        str2double(get(hObject,'String')) returns contents of EdANMthreshold as a double


% --- Executes during object creation, after setting all properties.
function EdAmax_max_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EdANMthreshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function EdCorrection_Callback(hObject, eventdata, handles)
% hObject    handle to EdCorrection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.NRthreshold = str2double(get(handles.EdCorrection,'String'));
handles.sys = set(handles.sys,'NRthreshold', handles.NRthreshold);
Residu_tot = handles.sys.R(handles.sys,handles.U0now);
if (norm(Residu_tot)>handles.NRthreshold)
    handles.U0now = NRcorrection(handles.sys, handles.U0now);
    handles.diagram.CurPoint.U0now=handles.U0now;
end
% %handles.diagram.Utnow=tangentvector(handles.sys,Jacobian(handles.sys,handles.U0now));
handles.Utnow=tangentvector(handles.sys,Jacobian(handles.sys,handles.U0now));
handles.hcurpoint=CurPointErase(handles.diagram,handles.params,handles.hcurpoint);
handles.hcurpoint=CurPointDisp(handles.diagram,handles.params,handles.sys);

fprintf('NR threshold is updated to %i ', handles.NRthreshold);
guidata(hObject,handles);

% Hints: get(hObject,'String') returns contents of EdCorrection as text
%        str2double(get(hObject,'String')) returns contents of EdCorrection as a double


% --- Executes during object creation, after setting all properties.
function EdCorrection_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EdCorrection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function EdCancelLast_Callback(hObject, eventdata, handles)
% hObject    handle to Ed_Hvalue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.CancelLast = str2double(get(handles.EdCancelLast,'String'));

guidata(hObject,handles);
% Hints: get(hObject,'String') returns contents of Ed_Hvalue as text
%        str2double(get(hObject,'String')) returns contents of Ed_Hvalue as a double


% --- Executes during object creation, after setting all properties.
function EdCancelLast_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Ed_Hvalue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in CkbPointdisplay.
function EdCorrType_Callback(hObject, eventdata, handles)
% hObject    handle to CkbCorrection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.NRmethod = str2double(get(handles.EdCorrType,'String'));
if (handles.NRmethod~=0) && (handles.NRmethod~=1) && (handles.NRmethod~=2)
    fprintf('The value for the type of correction must be 0 (off), 1 (fixed parameter) or 2 (variable parameter).');
    handles.NRmethod = 0;
end

handles.sys = set(handles.sys,'NRmethod', handles.NRmethod);

if handles.sys.NRmethod == 0
    fprintf('NR corrections off.');
elseif handles.sys.NRmethod == 1
    fprintf('NR corrections on, with fixed parameter.');
elseif handles.sys.NRmethod == 2
    fprintf('NR corrections on, with variable parameter.');
end

guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function EdCorrType_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EdCorrection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in Bifurcation.
function SimpleBifurcation_Callback(hObject, eventdata, handles)
% hObject    handle to Bifurcation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if get(handles.sys,'BifDetection') == 0
    handles.sys.BifDetection = 1;
    fprintf('Action. Bifurcation localization on');
else
    handles.sys.BifDetection = 0;
    fprintf('Action. Bifurcation localization off');
end

guidata(hObject, handles);
% Hint: get(hObject,'Value') returns toggle state of Bifurcation

% --- Executes on button press in Stability.
function Stability_Callback(hObject, eventdata, handles)
% hObject    handle to Stability (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if get(handles.sys,'StabilityCheck') == 0
    handles.sys.StabilityCheck = 1;
    fprintf('Action. Stability computation on');
else
    handles.sys.StabilityCheck = 0;
    fprintf('Action. Stability computation off');
end

guidata(hObject, handles);
% Hint: get(hObject,'Value') returns toggle state of Bifurcation


% --- Executes on button press in BtExit.
function BtExit_Callback(hObject, eventdata, handles)
% hObject    handle to BtExit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
confirm=  questdlg('Are you sure to close Manlab?','Close Manlab ?','Yes','No','No');
if strcmp(confirm,'Yes')
    for i=get(0,'Children')'
        if i~=handles.Manlab
            close(i);
        end
    end
    delete(handles.Manlab);
    clear all;
    clc
end


% --- Executes on button press in CkbPointdisplay.
function CkbPointdisplay_Callback(hObject, eventdata, handles)
% hObject    handle to CkbPointdisplay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if get(handles.CkbPointdisplay,'value')==1
    handles.params.pointdisplay = 1;
    fprintf('Action. point display on');
else
    handles.params.pointdisplay = 0;
    fprintf('Action. point display off');
end

guidata(hObject, handles);
% Hint: get(hObject,'Value') returns toggle state of CkbPointdisplay


% --- Executes on button press in CkbGlobaldisplay.
function CkbGlobaldisplay_Callback(hObject, eventdata, handles)
% hObject    handle to CkbGlobaldisplay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if get(handles.CkbGlobaldisplay,'value')==1
    handles.params.globaldisplay = 1;
    fprintf('Action. global display on');
else
    handles.params.globaldisplay = 0;
    fprintf('Action. global display off');
end

guidata(hObject, handles);
% Hint: get(hObject,'Value') returns toggle state of CkbGlobaldisplay

% --- Executes on button press in BtDisplaydiagram.
function BtDisplaydiagram_Callback(hObject, eventdata, handles)
% hObject    handle to BtDisplaydiagram (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
fprintf('Action: Global display of the diagram');

try
    handles.sys.global_display(handles.sys, handles.diagram.ChckPoint);
catch
    cellfun(@(sec)handles.sys.global_display(handles.sys,sec),handles.diagram.ChckPoint);
end
guidata(hObject, handles);
% --- Executes on button press in BtDisplaydiagram.


% --- Executes on mouse press over figure background.
function Manlab_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to Manlab (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


function EdJump_Callback(hObject, eventdata, handles)
% hObject    handle to EdJump (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.params.activecurve=str2double(get(handles.EdJump,'String'));

if (handles.params.activecurve <= length(handles.params.dispvars(:,1)))
    fprintf('Action: Active curve updated:  %i ',handles.params.activecurve);
else
    fprintf('Action: Provided number greater than the total number of curves  %i ',handles.params.activecurve);
    handles.params.activecurve=1;
end
guidata(hObject, handles);
% Hints: get(hObject,'String') returns contents of EdJump as text
%        str2double(get(hObject,'String')) returns contents of EdJump as a double


% --- Executes during object creation, after setting all properties.
function EdJump_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EdJump (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function ExportPoint_Callback(hObject, eventdata, handles)
% hObject    handle to ExportPoint (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global U
obj = CurPointSet(handles.diagram,handles.params,handles.sys);
U = obj.CurPoint.U0now;
disp(['Point exported in global variable U']);

function ExportSection_Callback(hObject, eventdata, handles)
% hObject    handle to ExportPoint (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global Section
Section = CheckPointExport(handles.diagram,'one');
disp(['Section exported in global variable Section']);

function ExportDiagram_Callback(hObject, eventdata, handles)
% hObject    handle to ExportPoint (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global Diagram
Diagram = CheckPointExport(handles.diagram,'all');
disp(['Diagram exported in global variable Diagram']);

function EdPerturbationSize_Callback(hObject, eventdata, handles)
% hObject    handle to EdPerturbationSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of EdPerturbationSize as text
%        str2double(get(hObject,'String')) returns contents of EdPerturbationSize as a double
handles.PerturbationSize = str2double(get(handles.EdPerturbationSize,'String'));
handles.sys = set(handles.sys,'PerturbationSize', handles.PerturbationSize);

handles.U0now = NRcorrection(handles.sys, get(handles.diagram,'U0now'));
handles.diagram.CurPoint.U0now=handles.U0now;

handles.Utnow=tangentvector(handles.sys,Jacobian(handles.sys,handles.U0now));
handles.hcurpoint=CurPointErase(handles.diagram,handles.params,handles.hcurpoint);
handles.hcurpoint=CurPointDisp(handles.diagram,handles.params,handles.sys);

if handles.PerturbationSize == 0
    handles.params.drawtype = {'-',':'};
    disp('The system is not perturbed.');
else
    handles.params.drawtype = {'--','-.'};
    fprintf('The system is perturbed with a perturbation vector proportional to %i ', handles.PerturbationSize);
end

guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function EdPerturbationSize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EdPerturbationSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function Manlab_SizeChangedFcn(hObject, eventdata, handles)
% hObject    handle to ExportPoint (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
