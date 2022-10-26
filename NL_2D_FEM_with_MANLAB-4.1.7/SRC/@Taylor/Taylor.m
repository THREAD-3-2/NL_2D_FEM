classdef Taylor
    %TAYLOR Taylor class constructor 
      %   t = Taylor(order,v) creates a Taylor object from the vector v
      %   containing the Taylor coefficients v0, {coef(j)}_{j=1,..,order}
      %   v0 is the value, coef(j) is the j^{th} Taylor coeffficient 
    properties
        value=0;
        order=0;
        coef=0;
    end
    
    %Class constructor
    methods
    
        function p=Taylor(varargin)

               if nargin==0
                    % EMPTY/DEFAULT CONSTRUCTOR
                    p.value=[];
                    p.coef=[];
                    p.order=[];
                    return;

                elseif nargin==1
                    obj=varargin{1};
                        switch class(obj)
                            case 'Taylor'
                                % COPY CONSTRUCTOR
                                    p.value=obj.value;
                                    p.coef=obj.coef;
                                    p.order=obj.order;
                                    return;
                            case 'double'
                                % CONVERT A 2-D or 1-D ARRAY TO TAYLORSERIES1
                                    p.value=obj;
                                    p.order=20;  % Default differentiation order
                                    p.coef=zeros([size(obj) p.order]);
                                    return;
                        end

                    elseif nargin==2                        
                        obj_order=varargin{1};
                        obj_value=varargin{2};
                        p.value=obj_value;
                        p.order=obj_order;
                        p.coef=zeros([size(obj_value) p.order]);
                        return;
                        
               elseif nargin==3
                        obj_order=varargin{1};
                        obj_value=varargin{2};
                        obj_coef=varargin{3};

                    % Checking if the dimensions of coef and value match well
                    scoef=size(obj_coef);
                svalue=size(obj_value);

                 if numel(scoef)==2|| scoef(1)~=svalue(1)...
                        ||scoef(2)~=svalue(2)|| scoef(3)~=obj_order
                        error('Objet.value and Objet.coef dimenssions mismatch')
                 end

                    p.value=obj_value;
                    p.order=obj_order;
                    p.coef=obj_coef;  
                    return;
               end
                error('Unsupported use of the Class TAYLORSERIES1.');
        end
    end
end
