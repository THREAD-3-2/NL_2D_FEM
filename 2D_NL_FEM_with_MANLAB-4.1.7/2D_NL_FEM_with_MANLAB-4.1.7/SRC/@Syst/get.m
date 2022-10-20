function [val] = get(sys,propName)
% @SYS/GET: Get sys properties from the specified object
% and return the value
switch propName
    case 'type'
        val = sys.type;
    case 'name'
        val = sys.name;
    case 'ninc'
        val = sys.ninc;
    case 'neq_tot'
        val = sys.neq_tot;
    case 'neq'
        val = sys.neq;
    case 'neq_aux'
        val = sys.neq_aux;
    case 'order'
        val = sys.order;
    case 'ANMthreshold'
        val = sys.ANMthreshold;
    case 'NRitemax'
        val = sys.NRitemax;
    case 'NRthreshold'
        val = sys.NRthreshold;
    case 'NRmethod'
        val = sys.NRmethod;
    case 'BifDetection'
        val = sys.BifDetection;
    case 'StabilityCheck'
        val = sys.StabilityCheck;
    case 'PerturbationSize'
        val = sys.PerturbationSize;
    case 'StabTol'
        val = sys.StabTol;
    case 'Amax_max'
        val = sys.Amax_max;
    case 'strinfos'
        val = sys.strinfos;
        
end


