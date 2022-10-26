function [F_full] = assemble_constant_vector(obj, type, varargin)
% assemble a vector depending on the provided type
% type: 'static_force_vector', 'periodic_force_vector', 'internal_force_vector', 'internal_force_vector_from_aux_var
% varargin depending on type:
%   static_force_vector: varargin empty
%   periodic_force_vector: varargin empty
%   internal_force_vector: varargin{1} = qs (given configuration)
%   internal_force_vector_from_aux_var: varargin{1} = aux_var (for MANLAB equation)

% initialize vector
if strcmp(type,'static_force_vector')
    F_full = zeros(3*obj.mesh.number_nodes, 1);
elseif strcmp(type,'periodic_force_vector')
    % find max number of forcing harmonic
    Hmax = 1;
    if isfield(obj.loads.periodic,'ponctual')
        for kk = 1:length(obj.loads.periodic.ponctual.node_list)
            harm = obj.loads.periodic.ponctual.node_list{kk}.harmonic;
            Hmax = max(Hmax, max(harm));
        end
    end
    if isfield(obj.loads.periodic,'distributed')
        harm = obj.loads.periodic.distributed.harmonic;
        if ~isempty(harm)
            Hmax = max(Hmax, max(harm));
        end
    end
    %  F_full is a concatenation of force vector, one for each harmonic   [F_1, F_2, ..]
    F_full = zeros(3*obj.mesh.number_nodes, Hmax);
elseif strcmp(type,'internal_force_vector')
    q_full = varargin{1}; % full vector of dof
    F_full = zeros(3*obj.mesh.number_nodes, size(q_full,2));
elseif strcmp(type,'internal_force_vector_from_aux_var')
    aux_var = varargin{1}; % full vector of dof
    F_full = zeros(3*obj.mesh.number_nodes, size(aux_var,2));
else
    error('Unknown vector type to assemble')
end
% Loop over the elements for distrubited forced of internal forces vectors
number_elements = obj.mesh.number_elements;
for i = 1:number_elements
    % Element node numbers
    nodeA  =  obj.mesh.connect(i,2);
    nodeB  =  obj.mesh.connect(i,3);
    % index of nodes dof
    index_global_A = 3*(nodeA - 1) + [1:3];
    index_global_B = 3*(nodeB - 1) + [1:3];
    index = [index_global_A index_global_B];
    % element information (length and angle w.r.t horizontal)
    [L_element, theta_element] = obj.elementOrientation(i);
    % rotation matrix to go from local frame to global frame
    rot_matrix = obj.transformationMatrix(theta_element);
    % elementary matrices computation
    if strcmp(type,'static_force_vector') % assemble disctributed force vector
        F_elem = zeros(6,1);
        if  isfield(obj.loads.static,'distributed')
            amp = obj.loads.static.distributed.amplitude;
            dir = obj.loads.static.distributed.direction;
            dir2 = [dir dir+3]; % dof for the 2 nodes of the element
            amp2 = [amp amp];   % amplitude for the 2 nodes of the element
            F_elem(dir2) = amp2*L_element/2; % set non zero components
        end
    elseif strcmp(type,'periodic_force_vector')
        
        if  isfield(obj.loads.periodic,'distributed')
            
            amp = obj.loads.periodic.distributed.amplitude;
            dir = obj.loads.periodic.distributed.direction;
            harm = obj.loads.periodic.distributed.harmonic;
            
            dir2 = [dir dir+3]; % dof for the 2 nodes of the element
            amp2 = [amp amp];   % amplitude for the 2 nodes of the element
            harm2 = [harm harm];   % amplitude for the 2 nodes of the element
            F_elem = zeros(6, max([1, harm]));
            for ii=1:length(dir2)
                F_elem(dir2(ii),harm2(ii)) = amp2(ii)*L_element/2; % set non zero components
            end
        else
            F_elem = zeros(6,1);
        end
    elseif strcmp(type,'internal_force_vector')
        qe = rot_matrix'*q_full(index,:);
        F_elem = rot_matrix*obj.elementary_internal_force_vector(qe, L_element); % rotate to go to global frame
    elseif strcmp(type,'internal_force_vector_from_aux_var')
        Fx = aux_var( 8*number_elements + i,:);
        Fy = aux_var( 9*number_elements + i,:) ;
        M = aux_var(10*number_elements + i,:) ;
        T2 = aux_var(11*number_elements + i,:) ;
        F_elem = rot_matrix*[-Fx;   % rotate to go to global frame
            -Fy;
            -M + T2*L_element/2;
            Fx;
            Fy;
            M + T2*L_element/2];
    end
    % assemble the element vector inside the full vector
    F_full(index,:) = F_full(index,:) + F_elem;
end
% Add point static and periodic forces
% do a loop on static ponctual forces:
if strcmp(type,'static_force_vector') % assemble point force vector
    if ~isempty(obj.loads.static.ponctual)
        if ~isempty(obj.loads.static.ponctual.node_list)
            for kk = 1:length(obj.loads.static.ponctual.node_list)
                if ~isempty(obj.loads.static.ponctual.node_list{kk})
                    % gather ppties of each forced node
                    node = obj.loads.static.ponctual.node_list{kk}.node;
                    dof = obj.loads.static.ponctual.node_list{kk}.dof;
                    amp = obj.loads.static.ponctual.node_list{kk}.amplitude;
                    % index of forced dof in full vector of dof
                    index = 3*(node-1)+dof;
                    F_full(index) = F_full(index) + amp';
                end
            end
        end
    end
end
% do a loop on periodic ponctual forces:
if strcmp(type,'periodic_force_vector') % assemble point force vector
    if ~isempty(obj.loads.periodic.ponctual)
        if ~isempty(obj.loads.periodic.ponctual.node_list)
            for kk = 1:length(obj.loads.periodic.ponctual.node_list)
                % gather ppties of each forced node
                if ~isempty(obj.loads.periodic.ponctual.node_list{kk})
                    node = obj.loads.periodic.ponctual.node_list{kk}.node;
                    dof = obj.loads.periodic.ponctual.node_list{kk}.dof;
                    amp = obj.loads.periodic.ponctual.node_list{kk}.amplitude;
                    harm = obj.loads.periodic.ponctual.node_list{kk}.harmonic;
                    % index of forced dof in full vector of dof
                    index = 3*(node-1)+dof;
                    % update force vector
                    for ii=1:length(index)
                        F_full(index(ii), harm(ii)) = F_full(index(ii), harm(ii)) + amp(ii);
                    end
                end
            end
        end
    end
end
end