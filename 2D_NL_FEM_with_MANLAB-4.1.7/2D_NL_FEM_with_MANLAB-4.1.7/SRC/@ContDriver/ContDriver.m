classdef ContDriver
    
    properties
   
        CurPoint   % a structure containing the Current Point of continuation  
                   % CurPoint is define in Manlab.m
                   % properties of CurPoint are : U0now,Utnow,Ut2,Status
                                 
                                
        ChckPoint = cell(0);    % a list of checkpoint

    end
    methods
          
         function obj = ContDriver(CurPoint)  
         
            obj.CurPoint = CurPoint;
            obj.ChckPoint = cell(0);  
            
        end
    end
end

