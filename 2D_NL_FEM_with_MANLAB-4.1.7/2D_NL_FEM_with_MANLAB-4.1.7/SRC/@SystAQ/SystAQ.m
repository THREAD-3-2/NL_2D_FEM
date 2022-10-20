classdef SystAQ < Syst
    
    properties
        
        iC=[] ; vC=[] ;                   % list for the sparse (order 1) tensor C
        iL=[] ; jL=[] ; vL=[] ;           % list for the sparse (order 2) tensor L
        iQ=[] ; jQ=[] ; kQ=[] ; vQ=[] ;   % list for the sparse (order 3) tensor Q
        idL=[]; jdL=[]; vdL=[];           % list for the sparse (order 2) tensor dL
        idQ=[]; jdQ=[]; kdQ=[]; vdQ=[];   % list for the sparse (order 3) tensor d
        
        ifunct=[] ; jfunct=[] ; kfunct=[] ; funct={}; % list for functions 
        %  equation 'ifunct'  is  U( jfunct) = funct ( U(kfunct) )
        
        writing;
       
    end
    
    methods
        
        function sys = SystAQ(nz,nz_aux,equations,point_display,global_display,parameters,writing,operators)
            
            if nargin<7; writing = 'standard'; end
            
            neq = nz;
            neq_aux = nz_aux;
            sys=sys@Syst('neq',neq,'neq_aux',neq_aux);
            
            sys.type = 'AQ';
            sys.writing = writing;
            
            sys.parameters = parameters;
            sys.equations = equations;
            sys.point_display = point_display;
            sys.global_display = global_display;
            
            if nargin<8
                switch writing
                    case 'standard'
                        %sys = get_operators_old(sys)
                        sys = get_operators(sys,sys.equations,sys.neq_tot,sys.ninc);
             
                    case 'vectorial'
                        disp('vectorial initialization of the system.');
                        %sys = get_operators_vec_old(sys)
                        sys = get_operators_vec(sys,sys.equations,sys.neq_tot,sys.ninc);
                    case 'subsystems'
                        disp('subsystem initialisation');
                    case 'vectorialFE'
                        disp('vectorial initialization of the system for FE model.');
                        sys = get_operators_vec_FE(sys);
                    case 'vectorialFE2'
                        disp('vectorial initialization of the system for FE2 model.');
                        sys = get_operators_vec_FE2(sys);
                end
            else
                sys.iC = operators.iC;   sys.vC = operators.vC;
                sys.iL = operators.iL;   sys.jL = operators.jL;   sys.vL = operators.vL;
                sys.iQ = operators.iQ;   sys.jQ = operators.jQ;   sys.kQ = operators.kQ; sys.vQ = operators.vQ;
                sys.idL = operators.idL; sys.jdL = operators.jdL; sys.vdL = operators.vdL;
                sys.idQ = operators.idQ; sys.jdQ = operators.jdQ; sys.kdQ = operators.kdQ; sys.vdQ = operators.vdQ;
                sys.ifunct = operators.ifunct; sys.jfunct = operators.jfunct; sys.kfunct = operators.kfunct;
                sys.funct  = operators.funct;
            end
            
            % affichage des listes operators
            %disp_operators(sys)
            
            
            % Arclength only on the main variables
            sys.arclengthdef = sparse((1:neq+1),ones(1,neq+1),1,sys.ninc,1);
            
            sys.R = @R;
            
            function [Rf] = R(sys,Uf)
                
                % C, L Q tensor
                Rf = sparse(sys.iC,ones(size(sys.iC)),sys.vC,sys.neq_tot,1) + ...
                     sparse(sys.iL,ones(size(sys.iL)),sys.vL.*Uf(sys.jL),sys.neq_tot,1) + ...
                     sparse(sys.iQ,ones(size(sys.iQ)),sys.vQ.*(Uf(sys.jQ).*Uf(sys.kQ)),sys.neq_tot,1);
                 
                for i=1:size(sys.ifunct,1)
                  Rf(sys.ifunct(i))= Uf(sys.kfunct(i))-sys.funct{i}(Uf(sys.jfunct(i)));
                end
                
            end
            
            function [ ] =disp_operators(sys)
                
              disp('Operators : C')
              if ~isempty(sys.iC)   
                for i=1:size(sys.iC)
                 disp([num2str(i),'  iC=',num2str(sys.iC(i)), ' vC=',num2str(sys.vC(i))])
                end
              end
              
              disp('Operators : L')
              if ~isempty(sys.iL)   
                for i=1:size(sys.iL)
                 disp([num2str(i),'  iL=',num2str(sys.iL(i)), '  jL=',num2str(sys.jL(i)),' vL=',num2str(sys.vL(i))])
                end
              end
              
              disp('Operators : Q')
              if ~isempty(sys.iQ)   
                for i=1:size(sys.iQ)
                 disp([num2str(i),'  iQ=',num2str(sys.iQ(i)), '  jQ=',num2str(sys.jQ(i)),'  kQ=',num2str(sys.kQ(i)),' vQ=',num2str(sys.vQ(i))])
                end
              end
              
              disp('Operators : dL')
              if ~isempty(sys.idL)   
                for i=1:size(sys.idL)
                 disp([num2str(i),'  idL=',num2str(sys.idL(i)), '  jdL=',num2str(sys.jdL(i)),' vdL=',num2str(sys.vdL(i))])
                end
              end
              
              disp('Operators : dQ')
              if ~isempty(sys.idQ)   
                for i=1:size(sys.idQ)
                 disp([num2str(i),'  idQ=',num2str(sys.idQ(i)), '  jdQ=',num2str(sys.jdQ(i)),'  kdQ=',num2str(sys.kdQ(i)),' vdQ=',num2str(sys.vdQ(i))]) 
                end 
              end
              
              disp('Operators : funct')
              if ~isempty(sys.ifunct)   
                for i=1:size(sys.ifunct)
                 disp([num2str(i),'  ifunct=',num2str(sys.ifunct(i)), '  jfunct=',num2str(sys.jfunct(i)),'  kfunct=',num2str(sys.kfunct(i))])
                end
              end
              
              
            end
        end
        
    end
end



