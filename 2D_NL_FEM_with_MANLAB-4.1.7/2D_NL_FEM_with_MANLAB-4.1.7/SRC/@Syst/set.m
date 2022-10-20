function sys = set(sys,varargin)
% @SYS/SET: Set sys properties and return the updated object
propertyArgin = varargin;
while length(propertyArgin) >=2
    prop = propertyArgin{1};
    val = propertyArgin{2};
    propertyArgin=propertyArgin(3:end);
    switch prop
        case 'name'
            sys.name = val;
            disp(['Problem: ' val]);
        case 'order'
            sys.order = val;
            disp(['series''s order = ' num2str(sys.order)]);
        case 'ANMthreshold'
            sys.ANMthreshold = val;
        case 'NRitemax'
            sys.NRitemax = val;
            disp(['Max number of iteration for NR corrections = ' num2str(sys.NRitemax)]);
        case 'NRthreshold'
            sys.NRthreshold = val;
        case 'NRmethod'
            sys.NRmethod = val;
        case 'BifDetection'
            sys.BifDetection = val;
            if val == 0; disp('Dectection of the Bifurcation off'); end
        case 'StabilityCheck'
            sys.StabilityCheck=val;
            if val == 1; disp('Stability computation on'); end
        case 'StabTol'
            sys.StabTol = val;
        case 'arclengthdef'
            if (length(val(1,:)) > 1)
                errordlg('Arclengthdef should be a column vector');
            end
            sys.arclengthdef = val;
        case 'PerturbationSize'
            sys.PerturbationSize = val;
        case 'type'
            sys.type = val;
        case 'ninc'
            sys.ninc = val;
        case 'neq_tot'
            sys.neq_tot = val;
        case 'neq'
            sys.neq = val;
        case 'neq_aux'
            sys.neq_aux = val;
        case 'Amax_max'
            sys.Amax_max = val;
    end
end


