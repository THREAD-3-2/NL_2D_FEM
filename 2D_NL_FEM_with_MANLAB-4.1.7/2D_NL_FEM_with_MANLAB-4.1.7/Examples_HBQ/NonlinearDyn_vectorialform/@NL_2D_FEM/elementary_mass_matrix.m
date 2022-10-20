function Me = elementary_mass_matrix(obj, L_element)
% elementary mass matrix of the timoshenko element
area =     obj.prop.area;
I =     obj.prop.inertia;
rho =  obj.prop.density;

% elementary mass matrix
Me =  L_element*rho/6*[ 2*area     0      0    area      0      0
                        0      2*area     0      0     area     0
                        0          0     2*I     0       0      I
                        area       0      0   2*area     0      0
                        0        area     0      0     2*area   0
                        0          0      I      0       0     2*I ];

end
