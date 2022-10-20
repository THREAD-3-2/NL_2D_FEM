geometricPoint = [    
                  1 -1 0 0;
                  2 1 0 0;
                  3 0 1 0 ;
                  4 0 1/2 0;
                  5 -1 -1 0;
                  6  1 -1 0;                  ];              
geometricElement = [1 1 2;
                    2 2 3; 
                    3 3 1 ;
                    4 3 4;
                    5 1 2 ; 
                    6 2 6;
                    7 6 5; 
                    8 5 1];
geometricDiscretisation = [5, 5, 5, 3, 5 ,5 5 5];

mesh = mesh_from_truss(geometricPoint, geometricElement, geometricDiscretisation);
fig = figure;
plot_deformed_mesh(mesh,zeros(mesh.number_nodes*3,1), fig, '-')
    
              