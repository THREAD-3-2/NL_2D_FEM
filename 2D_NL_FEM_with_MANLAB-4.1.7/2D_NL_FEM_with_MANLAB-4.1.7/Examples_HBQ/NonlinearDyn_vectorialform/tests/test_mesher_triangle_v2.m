geometricPoint = [1 -1 0 0;
                  2 1 0 0;
                  3 0 1 0 ;
                  4 0 1/2 0];              

geometricElement = [1 1 2;
                    2 2 3;
                    3 3 1 ;
                    4 3 4];
geometricDiscretisation = [5, 5, 5, 3];

mesh = mesh_from_truss(geometricPoint, geometricElement, geometricDiscretisation);
fig = figure;
plot_deformed_mesh(mesh,zeros(mesh.number_nodes*3,1), fig, '-')
    
              