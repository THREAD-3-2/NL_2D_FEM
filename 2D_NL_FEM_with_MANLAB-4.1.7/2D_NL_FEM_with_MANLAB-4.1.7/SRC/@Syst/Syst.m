classdef Syst
    
    properties        

        ninc  = 0;            % unknown number
        name  = '';
        order = 20;           % default serie truncature order
        ANMthreshold =1e-6;   % default ANM theshold
        NRitemax = 10;        % default max number of N-R iterations
        NRthreshold= 2e-5;    % default threshold for N-R corrections
        NRmethod = 0;         % N-R corrections off
        BifDetection = 1;     % Bifurcation detections on
        StabilityCheck = 0;   % Stability computation off
        StabTol = 1e-3;       % default tolerance for the stability 
        PerturbationSize = 0; % (global perturbation of the system) 's size
        Amax_max = 1e4;       % Maximum radius of convergence choosen.
        
        neq = 0;            % number of equations
        neq_aux = 0;        % number of auxiliar equations
        neq_tot = 0;        % total number of equations
        arclengthdef = 0;   % default path vector
        randvect=0;         % random vector used in the determination of the tangent vector
        pertvect=0;         % random positive vector used for computation of perturbed branches
        
        strinfos='';
        type=''
        
        equations;
        parameters;
        point_display;
        global_display;
        
        R;  % residu calculation (property for automatic differentiation)
    end
    
    methods
        function sys = Syst(varargin)
            % Constructor of the SYS class objects
            propertyArgin = varargin;
            while length(propertyArgin) >=2
                prop = propertyArgin{1};
                val=propertyArgin{2};
                propertyArgin=propertyArgin(3:end);
                switch prop
                    case 'neq'
                        sys.neq=val;
                    case 'neq_aux'
                        sys.neq_aux=val;
                    case 'name'
                        sys.name=val;
                    case 'order'
                        sys.order=val;
                    case 'ANMthreshold'
                        sys.ANMthreshold=val;
                    case 'NRthreshold'
                        sys.NRthreshold=val;
                    case 'NRmethod'
                        sys.NRmethod=val;
                    case 'NRitemax'
                        sys.NRitemax=val;
                    case 'BifDetection'
                        sys.BifDetection=val;
                    case 'StabilityCheck'
                        sys.StabilityCheck=val;
                    case 'StabTol'
                        sys.StabTol=val;
                    case 'type'
                        sys.type = val;
                    case 'parameters'
                        sys.parameters = val;
                    case 'point_display'
                        sys.point_display = val;
                    case 'global_display'
                        sys.global_display = val;
                    case 'R'
                        sys.R = val;
                end
            end
            
            sys.neq_tot  = sys.neq+sys.neq_aux;  % total number of equations
            sys.ninc = sys.neq_tot+1;                   % number of unknowns
            sys.arclengthdef = sparse(1:sys.neq+1,ones(1,sys.neq+1),ones(1,sys.neq+1),sys.ninc,1);    % default path vector
            sys.randvect = randn(1,sys.ninc);        % random vector used in the determination of the tangent vector for automatic differentiation only
            
            sys.pertvect = [rand(sys.neq,1);zeros(sys.neq_aux,1)];   % perturbation vector
            sys.pertvect = sys.pertvect/norm(sys.pertvect);          % normalization
            
            sys.strinfos = ['SYS : nequation=' num2str(sys.neq) ' nunknowns=' num2str(sys.ninc) ];
            disp(sys.strinfos);
        end
    end
end



