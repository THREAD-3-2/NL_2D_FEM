
classdef CheckPoint
    % These objects are created along the continuation process,
    % each time a new series is computed and stored.
    % The main property is the series U0(a=0) stored in U0 and its range of utility [0, Amax].
    % Other properties may be deduced from these two main properties
    properties
        
        % for storing the section
        U0   =0;    % the start point + series      (taylor type)
        Ut   =0;    % traveling direction at start point   (vector)
        Amax =0;    % range of utility of series
        Uend =0;    % the end point of the section  (vector)
        Utend=0;    % traveling direction at end      (vector)
        Upp  =0;    % U0 UI UII  ...  Uend, point by point representation of the section
        %  Upp(ninc,nbpts)  nbpts: number of point per section
        
        
        BifStatus={'nothing'};  % 'nothing' or 'simplebif' or 'Hopf' or 'Pitchfork' or 'NS' or 'PD' or 'Fold',...
        % for storing geometrical series and bifurcation data (when found)
        Ubif=0;       % bifurcation point
        Utf =0;       % tangent to the fondamental path
        Utb =0;       % tangent to the bifurcated path
        alpha=0;      % common ratio of the emerging geometrical series
        Uscale=0;     % scale vector of emerging geometrical series
        
        % for storing change of stability informations
        StabStatus={'stable','stable'};  % 'stable' or 'unstable' | beginning and end of the branch.
        ind_change;   % index of the point representation where the solution changes stability
        Ustab=0;      % change of stability point
        Astab=0;      % change of stability pseudo arc-length parameter
        Eigen;        % Stability analysis eigen values and eigen vectors
        Eigen_init;   % Eigen values/vectors at the beginning of the branch
        Eigen_end;    % Eigen values/vectors at the end of the branch
        drawtype_init='-';
        drawtype_end='-';
        
        % for displaying the projected bifurcation diagram
        dispcolors='';
        dispvars=0;
        %drawtype={'-',':'};
        %drawtype_pert = {'--','-.'};
        markers=''
        nbpts=0;
        ncurve=0;
        X=[];     % X coordinate of the point to be plotted on Proj Bif Diag
        Y=[];     % Y coordinate of the point to be plotted on Proj Bif Diag
        A=0;     % values a of the path variable of the point to be plotted
        Xbif=[];  % X coordinate of the bifurcation point (if one)
        Ybif=[];  % Y coordinate of the bifurcation point (if one)
        Xstab=[];  % X coordinate of the change of stability point (if one)
        Ystab=[];  % Y coordinate of the change of stability point (if one)
        
    end
    
    methods
        
        function obj = CheckPoint(U0,Ut,Amax,params,BifData,StabData)
            % construct  CheckPoint object from the two main properties U0, Amax
            % and BifData
            global Ck
            
            % parameters
            obj.dispcolors  = params.dispcolors;
            obj.dispvars    = params.dispvars;
            obj.markers     = params.markers;
            obj.nbpts       = params.nbpts;
            obj.ncurve      = size(obj.dispvars,1);
            
            %  Section definition
            obj.U0    = U0;         % the start point + series      (taylor type)
            obj.Ut    = Ut;         % traveling direction at start    (vector)
            obj.Amax  = Amax;       % range of utility of series
            obj.A     = zeros(obj.nbpts,1);
            obj.Upp   = zeros(size(Ut,1),obj.nbpts) ;
            for p=1:obj.nbpts
                a = (p-1)/(obj.nbpts-1)*Amax;
                obj.A(p) = a;                           % store values a of the path variable
                obj.Upp(:,p) = evalseries(U0,a,Ck); %  Upp = [U0 UI UII  ...  Uend]
            end
            
            obj.Uend  = obj.Upp(:,obj.nbpts) ;             % the end point                  (vector)
            obj.Utend = evalderiv(U0,obj.Amax,Ck) ;    % traveling direction at end      (vector)
            
            % store data  stability
            switch StabData.status{1}
                case 'stable'
                    obj.drawtype_init   = params.drawtype{1};
                case 'unstable'
                    obj.drawtype_init   = params.drawtype{2};
            end
            switch StabData.status{2}
                case 'stable'
                    obj.drawtype_end    = params.drawtype{1};
                case 'unstable'
                    obj.drawtype_end    = params.drawtype{2};
            end
            
            obj.BifStatus{1} = StabData.Eigen.type;
            obj.StabStatus = StabData.status;
            obj.Eigen = StabData.Eigen;
            obj.Eigen_init = StabData.Eigen_init;
            obj.Eigen_end = StabData.Eigen_end;
            obj.Ustab = StabData.Uchange;
            obj.Astab = StabData.Achange;
            
            [~,obj.ind_change] = find( (obj.A > obj.Astab)',1);
            
            % store data  bifurcation
            switch BifData.status
                case 'simplebif'
                    if strcmp(obj.BifStatus{1},'nothing') == 1
                        obj.BifStatus{1}=BifData.status;
                    else
                        obj.BifStatus{2}=BifData.status;
                    end
                    obj.alpha    =BifData.alpha;
                    obj.Uscale   =BifData.Uscale;
                    
                    obj.Ubif     =BifData.Ubif;
                    obj.Utf      =BifData.Utf;
                    obj.Utb      =BifData.Utb;
                    obj.Xbif=obj.Ubif(obj.dispvars(:,1));
                    obj.Ybif=obj.Ubif(obj.dispvars(:,2));
            end
            
            % data for the projected bifurcation diagram
            obj.X = obj.Upp(obj.dispvars(:,1),:)';   % X coordinate of the point to be plotted on Proj Bif Diag
            obj.Y = obj.Upp(obj.dispvars(:,2),:)';   % Y coordinate of the point to be plotted on Proj Bif Diag
            
            if obj.Astab > 0
                obj.Xstab=obj.Ustab(obj.dispvars(:,1));
                obj.Ystab=obj.Ustab(obj.dispvars(:,2));
            end
            
        end
    end
end

