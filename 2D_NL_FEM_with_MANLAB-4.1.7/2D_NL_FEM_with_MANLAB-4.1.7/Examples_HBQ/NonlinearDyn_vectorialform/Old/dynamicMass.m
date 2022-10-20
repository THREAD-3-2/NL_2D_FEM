function [M_element] = dynamicMass(L_element,params)
% makeMass creates the element mass matrix

density = params.density;
area = params.area;
I = params.I;

%%%% elemental mass matrix (lengths of each element could be different)
    
M_element = (L_element*density/6)*[ 2*area      0      0    area     0       0                                 
                                      0      2*area    0      0     area     0           
                                      0         0     2*I     0       0      I        
                                     area       0      0   2*area     0      0
                                      0        area    0      0     2*area   0
                                      0         0      I      0       0     2*I ];
end