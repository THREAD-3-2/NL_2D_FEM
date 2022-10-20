classdef NL_2D_FEM
    
    properties % (SetAccess=private)
        % FEM model Attributes
        geom     % geometry structure
        prop     % section and material properties structure
        mesh     % meh structure
        boundary % boundary condition structure
        loads    % static and periodic loads structures
        visu     % visualisation nodes structure 
        % FEM Solver Attributes
        solver = struct('type', 'homemade', 'tol_step', 1e-6, 'tol_res', 1e-6, 'n_iter_max', 500) % default solver parameters
        % FEM Equation Attributes        
        matrices = struct('mass', [], 'damping', [], 'stiffness_at_origin', [], 'stiffness_at_qs', [])
        vectors = struct('static_forces', [], 'periodic_forces', [], 'static_solution',[])        
    end
    
    methods
        %% automatic initialisation from model data structure
        function obj = set_model_from_mds(obj, mds)
            if isempty(mds.mesh.nodes) % mesh is empty, use provided geometry
                % set geometry (optional if a mesh is directly provided)
                obj = obj.set_geom(mds.geom.geom_node, mds.geom.geom_element, mds.geom.discretisation);
                % create mesh
                obj = obj.auto_mesh(); % auto mesh if geometry is given
            else
                obj = set_mesh(obj, mds.mesh.nodes, mds.mesh.connect);
            end
            % set visualized nodes
            obj = obj.set_visu(mds.visu.visu_node_list);
            % set properties
            obj = obj.set_prop(mds.prop.S, mds.prop.I, mds.prop.rho, mds.prop.E, mds.prop.poisson, mds.prop.k, mds.prop.alpha); % set propoerties (must be the same for all elements)
            % set boundary condition
            obj = obj.set_boundary(mds.boundary.bc_node_list);
            % set loads (forcing)
            %     static loads
            obj = obj.set_static_loads('ponctual', mds.force.static.static_ponctual_force_node_list);
            obj = obj.set_static_loads('distributed', mds.force.static.static_distributed_force_amplitude, mds.force.static.static_distributed_force_direction);
            %     dynamic loads
            obj = obj.set_periodic_loads('ponctual', mds.force.periodic.periodic_ponctual_force_node_list);
            obj = obj.set_periodic_loads('distributed', mds.force.periodic.periodic_distributed_force_amplitude, mds.force.periodic.periodic_distributed_force_direction, mds.force.periodic.periodic_distributed_force_harmonic);
            % assemble mass matrix and force vector
            obj = obj.initialise_matrices_and_vector();
            
        end
        
        %% manual initialisation methods
        function obj = set_geom(obj, geom_node, geom_element, discretisation)
            obj.geom.geom_node = geom_node;
            obj.geom.geom_element = geom_element;
            obj.geom.discretisation = discretisation;
            obj.mesh = obj.auto_mesh(); % auto mesh from given geometry
        end
        % visualisation
        function obj = set_visu(obj, visu_node_list)
            obj.visu.node_list = visu_node_list;
            index = [];
            for ii=1:length(visu_node_list)
                node = visu_node_list{ii}.node;   % integer
                dof = visu_node_list{ii}.dof; % subset of  {1,2,3}: list of variable size [1 2 3] or [2] or [2 3].
                % update list of global index
                index = [index 3*(node-1)+dof];
            end
            obj.visu.dof = index; % index of all visualized dof
        end
        % Properties        
        function obj = set_prop(obj, S, I, rho, E, eta, k, alpha) % all the same for every beam
            obj.prop.area = S;
            obj.prop.inertia = I;
            obj.prop.density = rho;
            obj.prop.young_mod = E;
            obj.prop.poisson = eta;
            obj.prop.shear_coeff_k = k;
            obj.prop.alpha = alpha;            
            % compute shear modulus from E and poisson
            obj.prop.shear_mod = E/(2*(1 + eta));
        end        
        %%      % meshing methods
        % auto mesh for given geometry
        function obj = auto_mesh(obj)
            [obj.mesh, mother_geometric_line] = mesh_from_truss(obj, obj.geom.geom_node, obj.geom.geom_element, obj.geom.discretisation);
        end
        % set data using user provided meh
        function obj = set_mesh(obj, nodes, connect)
            obj.mesh.nodes = nodes ;
            obj.mesh.connect = connect ;
            obj.mesh.number_nodes = size(nodes,1) ;
            obj.mesh.number_elements = size(connect,1) ;
            % obj.mesh.mother_geometric_line = 1; % TODO
        end        
        %% boundary condition methods
        function obj = set_boundary(obj, bc_node_list)
            obj.boundary.bc_node_list = bc_node_list;
            % construct index vector based on node / dof of bc condition
            prescribed_dof = [];
            for ii=1:length(bc_node_list)
                node = bc_node_list{ii}.node;   % integer
                dof = bc_node_list{ii}.dof; % subset of  {1,2,3}: list of variable size [1 2 3] or [2] or [2 3].
                % update list of global index
                prescribed_dof = [prescribed_dof, 3*(node-1) + dof]; % ndof_per_node * (node_number-1) + dof_number
            end
            dof_list = [1:3*obj.mesh.number_nodes];
            obj.boundary.dof_list = dof_list; % list of all dof number in the model
            obj.boundary.prescribed_dof = prescribed_dof; % list of prescribed dof number from b.c.
            obj.boundary.active_dof = setdiff(dof_list,prescribed_dof); % list of active dof number in the model
        end
        
        %% Loads methods
        % static loads
        function obj = set_static_loads(obj, type, varargin)
            
            if ~isempty(varargin)
                if strcmp(type,'ponctual')
                    obj.loads.static.ponctual.node_list = varargin{1};
                elseif strcmp(type,'distributed')
                    obj.loads.static.distributed.amplitude = varargin{1}; % ex: 10 or [10 10]
                    obj.loads.static.distributed.direction = varargin{2}; % ex: 1  or [1  2]
                else
                    error('undefined force type')
                end
            else
                if strcmp(type,'ponctual')
                    fprintf(1, 'No static point loads in the model \n')
                elseif strcmp(type,'distributed')
                    fprintf(1, 'No static distributed loads in the model \n')
                end
                
            end
        end
        % periodic loads
        function obj = set_periodic_loads(obj, type, varargin)
            if ~isempty(varargin)
                if strcmp(type,'ponctual')
                    obj.loads.periodic.ponctual.node_list = varargin{1};
                elseif strcmp(type,'distributed')
                    obj.loads.periodic.distributed.amplitude = varargin{1}; % ex: 10 or [10 10]
                    obj.loads.periodic.distributed.direction = varargin{2}; % ex: 1  or [1  2]
                    obj.loads.periodic.distributed.harmonic = varargin{3};    % ex: 1  or [2  5]
                else
                    error('undefined load type')
                end
            else
                if strcmp(type,'ponctual')
                    fprintf(1, 'No periodic point loads in the model \n')
                elseif strcmp(type,'distributed')
                    fprintf(1, 'No periodic distributed loads in the model \n')
                end
            end
        end
        % for updating static loads
        function obj = update_static_loads(obj, type, varargin)
            if ~isempty(varargin)
                if strcmp(type,'ponctual')
                    obj.loads.static.ponctual.node_list = varargin{1};
                elseif strcmp(type,'distributed')
                    obj.loads.static.distributed.amplitude = varargin{1}; % ex: 10 or [10 10]
                    obj.loads.static.distributed.direction = varargin{2}; % ex: 1  or [1  2]
                else
                    error('undefined force type')
                end
            else
                warning('no static force to update')
            end
            % reassemble the static force vector and update the model
            F_static_full  = assemble_constant_vector(obj, 'static_force_vector');
            obj.vectors.static_forces = F_static_full;
        end
        %% Solver setting
        function obj = set_solver(obj, solver)
                  obj.solver = solver;
        end        
        %% FEM Initialisation methods
        function obj = initialise_matrices_and_vector(obj)
            M_full = assemble_constant_matrix(obj, 'mass');
            K_full = assemble_constant_matrix(obj, 'stiffness_at_origin');
            F_static_full  = assemble_constant_vector(obj, 'static_force_vector');
            F_periodic_full  = assemble_constant_vector(obj, 'periodic_force_vector');
            % Matrices
            obj.matrices.mass = M_full;
            obj.matrices.damping = obj.prop.alpha*M_full;
            obj.matrices.stiffness_at_origin = K_full;
            % Vectors
            obj.vectors.static_forces = F_static_full;
            obj.vectors.periodic_forces = F_periodic_full;
            obj.vectors.null_vector = zeros(3*obj.mesh.number_nodes,1);
        end
    end
end

