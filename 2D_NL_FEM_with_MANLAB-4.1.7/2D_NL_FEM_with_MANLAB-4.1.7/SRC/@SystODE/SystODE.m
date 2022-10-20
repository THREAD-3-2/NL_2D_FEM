classdef SystODE < Syst
    
    properties
        
        H;  % Order of truncation of Fourier series
        
        subtype = 'autonomous' ; % subtype of system: {['autonomous'],'forced'}
        
        nz_tot ;    % number of equation of the DAE
        nz ;        % number of main equations
        nz_aux=0;   % number of auxiliary equations
        
        angfreq='omega';  % for forced system, the continuation parameter is omega by default.
        % Otherwise, sys.angfreq is the value of the
        % angfreq of the forcing terms
        
        % Small operator on the Algebro-differential equations
        id =[]; jd =[]; vd= [];                % list for the sparse (order 2) tensor d
        id1=[]; jd1=[]; vd1=[];                % list for the sparse (order 2) tensor d1
        idd=[]; jdd=[]; vdd=[];                % list for the sparse (order 2) tensor dd
        ic0=[]; vc0=[];                        % list for the sparse (order 1) tensor c0
        ic1=[]; vc1=[];                        % list for the sparse (order 1) tensor c1
        ic2=[]; vc2=[];                        % list for the sparse (order 1) tensor c2
        il0=[]; jl0=[]; vl0=[];                % list for the sparse (order 2) tensor l0
        il1=[]; jl1=[]; vl1=[];                % list for the sparse (order 2) tensor l1
        iq =[]; jq =[]; kq =[]; vq=[];         % list for the sparse (order 3) tensor q
        uiq=[];i_iq=[];
        
        iforce0=[]; hforce0=[]; vforce0=[];       % list for the sparse (order 2) tensor force
        iforce1=[]; hforce1=[]; vforce1=[];       % list for the sparse (order 2) tensor force
        iforce2=[]; hforce2=[]; vforce2=[];       % list for the sparse (order 2) tensor force
        
        idl=[]; jdl=[]; vdl=[];                % list for the sparse (order 2) tensor dl
        idq=[]; jdq=[]; kdq=[]; vdq=[];        % list for the sparse (order 3) tensor dq
        
        zi_phase1=1; % index of the state variable receiving the first phase condition.
        zi_phase2=2; % index of the state variable receiving the second phase condition.
        
        writing = 'standard';
    end
    
    methods
        
        function sys = SystODE(nz,nz_aux,H,equations,point_display,global_display,parameters,subtype,writing)
            
            if nargin<8; subtype = 'autonomous'; end
            if nargin<9; writing = 'standard'; end
            
            if sum(H) == 0
                neq = nz;
                neq_aux = nz_aux;
                type = 'Equ';
            elseif (numel(H) == 1 || H(1) == 0 || H(2) == 0)
                H = sum(H);
                neq = nz*(2*H+1)+1;
                neq_aux = nz_aux*(2*H+1)+2;
                type = 'HBM';
            else
                nb_coef = ((2*H(2)+1)*H(1) + H(2))*2 +1;
                if H(1) ~= H(2); errordlg('Different number of harmonics is not implemented yet.'); end
                neq = nz*nb_coef+2;
                neq_aux = nz_aux*nb_coef+5;
                type = 'QPHBM';
            end
            
            sys = sys@Syst('neq',neq,'neq_aux',neq_aux);
            
            sys.type = type;
            sys.subtype = subtype;
            sys.writing = writing;
            sys.R = @R;
            
            sys.parameters = parameters;
            sys.equations = equations;
            sys.point_display = point_display;
            sys.global_display = global_display;
            
            sys.nz = nz;
            sys.nz_aux = nz_aux;
            sys.nz_tot = nz + nz_aux;
            sys.H  = H(1) ;    % harmonic number
            
            % Arclength only on the main variables
            sys.arclengthdef = sparse(1:sys.neq+1,ones(sys.neq+1,1),1,sys.ninc,1);
            
            % angfreq of the oscillations for forced systems.
            if strcmp(sys.subtype,'forced')
                try sys.angfreq = sys.parameters.angfreq; % the angfreq of the forcing terms is constant.
                catch; sys.angfreq = 'omega'; % the angfreq of the forcing terms is the continuation parameter.
                end
            end
            
            % Creation of the little operators
            switch writing
                case 'vectorial'
                    sys = sys.get_operators_vec;
                case 'standard'
                    sys = sys.get_operators;
                otherwise
                    error("SystODE : Uncorrect writing type of the system. It should be either 'standard' or 'vectorial'.");
            end
            
            if strcmp(type,'QPHBM') && ~isempty(sys.idl)
                warndlg('SystODE : dL,dQ does not work for Quasi-periodic solutions in the general case.');
            end
            
            if ~isempty(find(sys.jd(sys.vd~=0)>sys.nz,1))
                warndlg("SystODE : Stability computation not available when time derivatives of auxiliary variables appear.");
                return
            end
            
            % Specific treatment of the products to avoid a long loop in Fpnl.m
            [sys.uiq,sys.i_iq,~] = unique(sys.iq); % Select one product per equation to allow a vectorial allocation
            
            
            %% Residue functions
            function [Rtot] = R(sys,Utot)
                switch sys.type
                    case 'Equ'
                        Rtot=R_Equ(sys,Utot);
                    case 'HBM'
                        Rtot=R_HBM(sys,Utot);
                    case 'QPHBM'
                        Rtot=R_QPHBM(sys,Utot);
                end
            end
            
            function [Rtot] = R_Equ(sys,Utot)
                %  Compute  Rtot(U)= C + L(U)   + Q(U,U)   + dL(dU) + Qh(U,dU);
                
                lambda = Utot(sys.nz+1);
                Uvar = Utot([1:sys.nz sys.nz+2:sys.nz_tot+1]);
                
                Ualg = [Uvar;lambda];
                
                Rtot = sys.equations(sys,0,Ualg,zeros(size(Ualg)),zeros(size(Ualg)));
                
                %                 Rtot =sparse(sys.ic0,ones(1,size(sys.ic0,1)),sys.vc0,sys.nz_tot,1) ...
                %                     + sparse(sys.ic1,ones(1,size(sys.ic1,1)),sys.vc1*lambda,sys.nz_tot,1) ...
                %                     + sparse(sys.ic2,ones(1,size(sys.ic2,1)),sys.vc2*lambda^2,sys.nz_tot,1) ...
                %                     + sparse(sys.il0,ones(1,size(sys.il0,1)),sys.vl0.*Uvar(sys.jl0),sys.nz_tot,1) ...
                %                     + sparse(sys.il1,ones(1,size(sys.il1,1)),(lambda*sys.vl1).*Uvar(sys.jl1),sys.nz_tot,1) ...
                %                     + sparse(sys.iq,ones(1,size(sys.iq,1)),sys.vq.*(Uvar(sys.jq).*Uvar(sys.kq)),sys.nz_tot,1) ...
                %                     + sparse(sys.idl,ones(1,size(sys.idl,1)),sys.vdl.*Uvar(sys.jdl),sys.nz_tot,1);
            end
            
            function [Rtot] = R_HBM(sys,Utot)
                %  Compute  R(U)
                H    = sys.H;   % Harmonic number
                nz_tot   = sys.nz_tot;  % number of equation of the DAE system
                DHp1 = 2*H+1 ;
                
                % Extract informations from Utot
                [Ztot,omega,lambda,omega2,lambdaomega] = sys.get_Ztot(Utot);
                Rfourier=zeros(DHp1,nz_tot);
                
                %   R = C0 + lambda C1  + L0(U) + lambda L1(U) +
                %       omega D(U) + lambdaomega D1(U) + omega2 DD(U) + Q(U,U)
                %
                %   R is formed by block using the list ic0, vc0, ......
                for i=1:numel(sys.iforce0)
                    Rfourier(sys.hforce0(i),sys.iforce0(i)) = Rfourier(sys.hforce0(i),sys.iforce0(i)) - sys.vforce0(i);
                end
                
                for i=1:numel(sys.iforce1)
                    Rfourier(sys.hforce1(i),sys.iforce1(i)) = Rfourier(sys.hforce1(i),sys.iforce1(i)) - lambda*sys.vforce1(i);
                end
                
                for i=1:numel(sys.iforce2)
                    Rfourier(sys.hforce2(i),sys.iforce2(i)) = Rfourier(sys.hforce2(i),sys.iforce2(i)) - lambda^2*sys.vforce2(i);
                end
                
                for i=1:numel(sys.ic0)
                    Rfourier(1,sys.ic0(i))= Rfourier(1,sys.ic0(i)) + sys.vc0(i) ;
                end
                
                for i=1:numel(sys.ic1)
                    Rfourier(1,sys.ic1(i))= Rfourier(1,sys.ic1(i)) + sys.vc1(i)*lambda ;
                end
                
                for i=1:numel(sys.ic2)
                    Rfourier(1,sys.ic2(i))= Rfourier(1,sys.ic2(i)) + sys.vc2(i)*lambda^2 ;
                end
                
                for i=1:numel(sys.il0)
                    Rfourier(:,sys.il0(i))= Rfourier(:,sys.il0(i)) + sys.vl0(i)*Ztot(:,sys.jl0(i)) ;
                end
                
                for i=1:numel(sys.il1)
                    Rfourier(:,sys.il1(i))= Rfourier(:,sys.il1(i)) + (sys.vl1(i)*lambda)*Ztot(:,sys.jl1(i)) ;
                end
                
                for i=1:numel(sys.id)
                    Rfourier(:,sys.id(i))= Rfourier(:,sys.id(i)) + (sys.vd(i)*omega)*sys.D(Ztot(:,sys.jd(i))) ;
                end
                
                for i=1:numel(sys.id1)
                    Rfourier(:,sys.id1(i))= Rfourier(:,sys.id1(i)) + (sys.vd1(i)*lambdaomega)*sys.D(Ztot(:,sys.jd1(i))) ;
                end
                
                for i=1:numel(sys.idd)
                    Rfourier(:,sys.idd(i))= Rfourier(:,sys.idd(i)) + (sys.vdd(i)*omega2)*sys.D(sys.D(Ztot(:,sys.jdd(i)))) ;
                end
                
                for i=1:numel(sys.iq)
                    Rfourier(:,sys.iq(i))= Rfourier(:,sys.iq(i)) + sys.vq(i)*( sys.Prod(Ztot(:,sys.jq(i)),Ztot(:,sys.kq(i))) ) ;
                end
                
                Ztot_t0 = sum(Ztot(1:H+1,:),1);
                
                % phase equation
                switch sys.subtype
                    case 'autonomous'
                        Rphase= (1:H)*Ztot(H+2:DHp1,sys.zi_phase1); % derivative of zi_phase1 variable null at t=0.
                        %Rphase= Ztot(H+2,sys.zi_phase1);            %first sine of zi_phase1 variable null.
                    case 'forced'
                        if isfloat(sys.angfreq)
                            Rphase = omega - sys.angfreq;
                        else
                            Rphase = omega - lambda;
                        end
                end
                
                Romega2 = omega2 - omega*omega;
                Rlambdaomega = lambdaomega - lambda*omega;
                
                if ~isempty(sys.idl)
                    for i=1:numel(sys.idl)
                        Rfourier(:,sys.idl(i))= Rfourier(:,sys.idl(i)) + sys.vdl(i)*sys.D(Ztot(:,sys.jdl(i))) ;
                    end
                    
                    for i=1:numel(sys.idq)
                        Rfourier(:,sys.idq(i))= Rfourier(:,sys.idq(i)) + sys.vdq(i)*( sys.Prod(Ztot(:,sys.jdq(i)),sys.D(Ztot(:,sys.kdq(i)))) ) ;
                    end
                    
                    cond_ini = sys.equations(sys,0,[Ztot_t0 lambda]',zeros(sys.nz_tot,1),zeros(sys.nz_tot,1));
                    Rfourier(1,sys.idl) = cond_ini(sys.idl);
                end
                
                % We divide Rfourier between the main equations and the auxiliary
                % equations. We add the phase equation and the definition of omega^2 and lambda*omega
                Rfourier = reshape(Rfourier,nz_tot*DHp1,1);
                R = [Rfourier(1:sys.neq-1); Rphase];
                
                Raux = [Rfourier(sys.neq:end);Romega2;Rlambdaomega];
                Rtot = [R;Raux];
                
            end
            
            function [Rtot] = R_QPHBM(sys,Utot)
                %  Compute  R(U)
                H    = sys.H;   % Harmonic number
                nz_tot   = sys.nz_tot;  % number of equation of the DAE system
                coef_per_var = (H+1)*4*H+1;
                
                % Extract informations from Utot
                [Ztot,omega1,omega2,lambda,omega1sq,omega2sq,omega1omega2,lambdaomega1,lambdaomega2] = get_Ztot(sys,Utot);
                    Ztot_t0 = sum(Ztot(1:(coef_per_var-1)/2+1,:),1);

                    Rfourier=zeros(coef_per_var,nz_tot);
                
                %   R = C0 + lambda C1  + L0(U) + lambda L1(U) +
                %       omega D(U) + lambdaomega D1(U) + omega2 DD(U) + Q(U,U)
                %
                %   R is formed by block using the list ic0, vc0, ......
                for i=1:size(sys.iforce0,1)
                    Rfourier( (sys.hforce0(i)-1)*(2*H+1)+1,sys.iforce0(i)) = Rfourier( (sys.hforce0(i)-1)*(2*H+1)+1,sys.iforce0(i)) - sys.vforce0(i);
                end
                
                for i=1:size(sys.iforce1,1)
                    Rfourier( (sys.hforce1(i)-1)*(2*H+1)+1,sys.iforce1(i)) = Rfourier( (sys.hforce1(i)-1)*(2*H+1)+1,sys.iforce1(i)) - lambda*sys.vforce1(i);
                end
                
                for i=1:size(sys.iforce2,1)
                    Rfourier( (sys.hforce2(i)-1)*(2*H+1)+1,sys.iforce2(i)) = Rfourier( (sys.hforce2(i)-1)*(2*H+1)+1,sys.iforce2(i)) - lambda^2*sys.vforce2(i);
                end
                
                for i=1:size(sys.ic0,1)
                    Rfourier(1,sys.ic0(i))= Rfourier(1,sys.ic0(i)) + sys.vc0(i) ;
                end
                
                for i=1:size(sys.ic1,1)
                    Rfourier(1,sys.ic1(i))= Rfourier(1,sys.ic1(i)) + sys.vc1(i)*lambda ;
                end
                
                for i=1:size(sys.ic2,1)
                    Rfourier(1,sys.ic2(i))= Rfourier(1,sys.ic2(i)) + sys.vc2(i)*lambda^2 ;
                end
                
                for i=1:size(sys.il0,1)
                    Rfourier(:,sys.il0(i))= Rfourier(:,sys.il0(i)) + sys.vl0(i)*Ztot(:,sys.jl0(i)) ;
                end
                
                for i=1:size(sys.il1,1)
                    Rfourier(:,sys.il1(i))= Rfourier(:,sys.il1(i)) + (sys.vl1(i)*lambda)*Ztot(:,sys.jl1(i)) ;
                end
                
                for i=1:size(sys.id,1)
                    Rfourier(:,sys.id(i))= Rfourier(:,sys.id(i)) + sys.vd(i)*sys.D(Ztot(:,sys.jd(i)),omega1,omega2) ;
                end
                
                for i=1:size(sys.id1,1)
                    Rfourier(:,sys.id1(i))= Rfourier(:,sys.id1(i)) + sys.vd1(i)*sys.D(Ztot(:,sys.jd1(i)),lambdaomega1,lambdaomega2) ;
                end
                
                for i=1:size(sys.idd,1)
                    Rfourier(:,sys.idd(i))= Rfourier(:,sys.idd(i)) + sys.vdd(i)*sys.DD(Ztot(:,sys.jdd(i)),omega1sq,omega2sq,omega1omega2) ;
                end
                
                for i=1:size(sys.iq,1)
                    Rfourier(:,sys.iq(i))= Rfourier(:,sys.iq(i)) + sys.vq(i)*( sys.Prod(Ztot(:,sys.jq(i)),Ztot(:,sys.kq(i))) ) ;
                end
                
                
                % phase equation
                switch sys.subtype
                    case 'autonomous'
                        Rphase1= Ztot((coef_per_var-1)/2+1+(1+2*H),sys.zi_phase1);        % sin(omega1) of variable "zi_phase1"
                        %Rphase1 = Ztot_t0(sys.zi_phase1);
                    case 'forced'
                        if isfloat(sys.angfreq)
                            Rphase1 = omega1 - sys.angfreq;
                        else
                            Rphase1 = omega1 - lambda;
                        end
                end
                Rphase2= Ztot((coef_per_var-1)/2+2,sys.zi_phase2); % sin(omega2) de la variable zi_phase2
                Rphase = [Rphase1;Rphase2];
                
                Romega1sq = omega1sq - omega1*omega1;
                Romega2sq = omega2sq - omega2*omega2;
                Romega1omega2 = omega1omega2 - omega1*omega2;
                Rlambdaomega1 = lambdaomega1 - lambda*omega1;
                Rlambdaomega2 = lambdaomega2 - lambda*omega2;
                
                
                if ~isempty(sys.idl)

                    % Assumption : the variable after kdq(i) is the derivative of the variable kdq(i)
                    for i=1:size(sys.idl,1)
                        Rfourier(:,sys.idl(i))= Rfourier(:,sys.idl(i)) + sys.vdl(i)*sys.D(Ztot(:,sys.jdl(i)),omega1,omega2) ;
                    end
                    % Assumption : the variable after kdq(i) is the derivative of the variable kdq(i)
                    for i=1:size(sys.idq,1) % Assumption : the variable after kdq(i) is the derivative of the variable kdq(i)
                        Rfourier(:,sys.idq(i))= Rfourier(:,sys.idq(i)) + sys.vdq(i)*sys.Prod(Ztot(:,sys.jdq(i)),Ztot(:,sys.kdq(i)+1));
                    end
                    
                    cond_ini = sys.equations(sys,0,[Ztot_t0 lambda]',zeros(sys.nz_tot,1),zeros(sys.nz_tot,1));
                    Rfourier(1,sys.idl) = cond_ini(sys.idl);
                end
                
                % We divide Rfourier between the main equations and the auxiliary
                % equations. We add the phase equation and the definition of omega^2 and lambda*omega
                Rfourier = reshape(Rfourier,nz_tot*coef_per_var,1);
                R = [Rfourier(1:sys.neq-2); Rphase];
                
                Raux = [Rfourier(sys.neq-1:end);Romega1sq;Romega2sq;Romega1omega2;Rlambdaomega1;Rlambdaomega2];
                Rtot = [R;Raux];
            end
            
        end
        
        %% Methods for initialization
        function [Unew] = init_Hdiff(sys,Uold)
            
            switch sys.type
                case 'HBM'
                    [Unew] = sys.init_Hdiff_HBM(Uold);
                case 'QPHBM'
                    [Unew] = sys.init_Hdiff_QPHBM(Uold);
            end
            
        end
        
        function [Unew] = init_Hdiff_HBM(sys,Uold)
            
            Hold = ((size(Uold,1)-4)/sys.nz_tot-1)/2;
            
            neq_old = sys.nz*(2*Hold+1);
            neq_aux_old = sys.nz_aux*(2*Hold+1);
            Zold = reshape(Uold([1:neq_old neq_old+3:neq_old+3+neq_aux_old-1]) , 2*Hold+1 , sys.nz_tot);
            omega = Uold(neq_old+1);
            lambda = Uold(neq_old+2);
            
            Hnew = sys.H;
            DHp1 = 2*Hnew+1;
            
            if Hold >= Hnew
                Znew = Zold([(1:Hnew+1) (Hold+2:Hold+1+Hnew)],:);
            else
                Znew = zeros(DHp1,sys.nz_tot);
                Znew([(1:Hold+1) (Hnew+2:Hnew+1+Hold)],:) = Zold;
            end
            
            Unew = sys.init_U0(Znew,omega,lambda);
            
        end
        
        function [Unew] = init_Hdiff_QPHBM(sys,Uold)
            
            Hold = floor(sqrt(((size(Uold,1)-8)/sys.nz_tot-1)/4));
            nb_coef_old = 2*((2*Hold+1)*Hold + Hold)+1;
            
            neq_old = sys.nz*nb_coef_old;
            neq_aux_old = sys.nz_aux*nb_coef_old;
            Zold = reshape(Uold([1:neq_old neq_old+4:neq_old+4+neq_aux_old-1]) , nb_coef_old , sys.nz_tot);
            omega1 = Uold(neq_old+1);
            omega2 = Uold(neq_old+2);
            lambda = Uold(neq_old+3);
            
            Hnew = sys.H;
            nb_coef_new = 2*((2*Hnew+1)*Hnew + Hnew)+1;
            
            if Hnew >= Hold
                Znew = zeros(nb_coef_new,sys.nz_tot);
                indices_cos = [(1:Hold+1)' ; reshape(repmat(2*Hnew+1 + 1 + (-Hold:1:Hold)',1,Hold) + (0:1:Hold-1)*(2*Hnew+1),(2*Hold+1)*Hold,1)];
                indices_sin = indices_cos(2:end)+(2*Hnew+1)*Hnew + Hnew;
                Znew([indices_cos;indices_sin],:) = Zold;
            else
                indices_cos = [(1:Hnew+1)' ; reshape(repmat(2*Hold+1 + 1 + (-Hnew:1:Hnew)',1,Hnew) + (0:1:Hnew-1)*(2*Hold+1),(2*Hnew+1)*Hnew,1)];
                indices_sin = indices_cos(2:end)+(2*Hold+1)*Hold + Hold;
                Znew = Zold([indices_cos;indices_sin],:);
            end
            
            Unew = sys.init_U0(Znew,omega1,omega2,lambda);
            
        end
        
        function [U0] = init_U0(sys,Z,varargin)
            if nargin == 4
                U0 = sys.init_U0_HBM(Z,varargin{:});
            elseif nargin == 5
                U0 = sys.init_U0_QPHBM(Z,varargin{:});
            else
                error('SystODE/init_U0.m : Wrong number of argument.');
            end
        end
        
        function [U0] = init_U0_HBM(sys,Z,omega,lambda)
            
            DHp1 = 2*sys.H+1;
            
            if numel(Z) == DHp1*sys.nz_tot
                Ufourier = reshape(Z(:,1:sys.nz),DHp1*sys.nz,1);
                Uaux_fourier = reshape(Z(:,sys.nz+1:end),DHp1*sys.nz_aux,1);
            end
            
            U0 = [Ufourier;omega;lambda;Uaux_fourier;omega^2;lambda*omega];
            
        end
        
        function [U0] = init_U0_QPHBM(sys,Z,omega1,omega2,lambda)
            
            nb_coef = 2*((2*sys.H+2)*sys.H)+1;
            
            Ufourier = reshape(Z(:,1:sys.nz),nb_coef*sys.nz,1);
            Uaux_fourier = reshape(Z(:,sys.nz+1:end),nb_coef*sys.nz_aux,1);
            
            U0 = [Ufourier;omega1;omega2;lambda;Uaux_fourier;omega1^2;omega2^2;omega1*omega2;lambda*omega1;lambda*omega2];
            
        end
        
        function [coord] = getcoord(sys,variablename,varargin)
            switch sys.type
                case 'HBM'
                    coord = getcoord_HBM(sys,variablename,varargin{:});
                case 'QPHBM'
                    coord = getcoord_QPHBM(sys,variablename,varargin{:});
                otherwise
                    error('@SystODE/getcoord.m : Undefined function for this type of system.');
            end
        end
        
        function [coord] = getcoord_HBM(sys,variablename,varargin)
            % function [coord] = getcoord(sys,variablename,varargin)
            % gives the coordinates of the variable 'variablename' in the total vector
            % of unknowns.
            %
            % The third and fourth argument values specifies the number of the variable
            % and of the harmonics considered.
            %
            % Example : getcoord('lambda') would give coord = sys.neq+1.
            %           getcoord('cos',1,5) would give the index of the coefficient =
            % of the fifth cosine of the first variable : that is 6.
            
            DHp1 = 2*sys.H+1;
            
            switch variablename
                case 'omega'
                    coord = sys.neq;
                case 'lambda'
                    coord = sys.neq+1;
                case 'cos'
                    ivar = varargin{1};
                    ih = varargin{2};
                    if ivar <= sys.nz
                        coord = (ivar-1)*DHp1 + ih+1;
                    else
                        coord = sys.neq + 1 + (ivar-sys.nz-1)*DHp1 + ih+1;
                    end
                case 'sin'
                    ivar = varargin{1};
                    ih = varargin{2};
                    if ivar <= sys.nz
                        coord = (ivar-1)*DHp1 + sys.H + ih+1;
                    else
                        coord = sys.neq + 1 + (ivar-sys.nz-1)*DHp1 + sys.H + ih+1;
                    end
            end
        end
        
        function [coord] = getcoord_QPHBM(sys,variablename,varargin)
            % function [coord] = getcoord(sys,variablename,varargin)
            % gives the coordinates of the variable 'variablename' in the total vector
            % of unknowns.
            %
            % The third and fourth argument values specifies the number of the variable
            % and of the harmonics considered.
            %
            % Example : getcoord('lambda') would give coord = sys.neq+1.
            %           getcoord('cos',1,[0 5]) would give the index of the coefficient =
            % of cos(5*omega2*t) of the first variable : that is 6.
            
            %DHp1 = 2*sys.H(2)+1;
            %coefpervar = 2*(sys.H(1)*DHp1+sys.H(2))+1;
            
            DHp1 = 2*sys.H+1;
            coefpervar = 2*(sys.H*DHp1+sys.H)+1;
            
            switch variablename
                case {'omega','omega1'}
                    coord = sys.neq-1;
                case 'omega2'
                    coord = sys.neq;
                case 'lambda'
                    coord = sys.neq+1;
                case 'cos'
                    ivar = varargin{1};
                    ih = varargin{2};
                    if ivar <= sys.nz
                        coord = (ivar-1)*coefpervar + ih(1)*DHp1 + ih(2)+1;
                    else
                        coord = sys.neq + 1 + (ivar-sys.nz-1)*coefpervar + ih(1)*DHp1 + ih(2)+1;
                    end
                case 'sin'
                    ivar = varargin{1};
                    ih = varargin{2};
                    if ivar <= sys.nz
                        coord = (ivar-1)*coefpervar + (coefpervar-1)/2 + ih(1)*DHp1 + ih(2)+1;
                    else
                        coord = sys.neq + 1 + (ivar-sys.nz-1)*coefpervar + (coefpervar-1)/2 + ih(1)*DHp1 + ih(2)+1;
                    end
            end
            
        end
        
        %% Methods specific to HBM or QPHBM types.
        function [Ztot,varargout] = get_Ztot(sys,Utot)
            
            switch sys.type
                case 'HBM'
                    DHp1 = sys.H*2+1;
                    Z      = reshape(Utot(1:sys.neq-1),DHp1,sys.nz);
                    omega  = Utot(sys.neq);
                    lambda = Utot(sys.neq+1) ;
                    
                    omegasq = Utot(end-1);
                    lambdaomega =Utot(end);
                    Zaux   = reshape(Utot(sys.neq+2:end-2),DHp1,sys.nz_aux);
                    
                    Ztot = [Z Zaux];
                    varargout = {omega,lambda,omegasq,lambdaomega};
                case 'QPHBM'
                    
                    coef_per_var = (sys.H+1)*4*sys.H+1;
                    Z      = reshape(Utot(1:sys.neq-2),coef_per_var,sys.nz);
                    omega1 =Utot(sys.neq-1);
                    omega2 =Utot(sys.neq);
                    lambda=Utot(sys.neq+1) ;
                    
                    omega1sq = Utot(end-4);
                    omega2sq = Utot(end-3);
                    omega1omega2 = Utot(end-2);
                    lambdaomega1 = Utot(end-1);
                    lambdaomega2 = Utot(end);
                    Zaux = reshape(Utot(sys.neq+2:end-5),coef_per_var,sys.nz_aux);
                    
                    Ztot = [Z , Zaux];
                    varargout = {omega1,omega2,lambda,omega1sq,omega2sq,omega1omega2,lambdaomega1,lambdaomega2};
            end
        end
        
        function DU = D(sys,U,varargin)
            switch sys.type
                case 'HBM'
                    %  Compute  D.U
                    Hu=(size(U,1)-1)/2   ;   % harmonic number
                    if nargin < 3
                        DU = [ zeros(1,size(U,2)) ; (1:Hu)'.*U(Hu+2:end,:) ; (-(1:Hu))'.*U(2:Hu+1,:) ];
                    else
                        DU = 1i*((-Hu:Hu)'.*U);
                    end
                case 'QPHBM'
                    %  Compute  D.U
                    omega1 = varargin{1};
                    omega2 = varargin{2};
                    Hu=floor(sqrt((size(U,1)-1)/4));   % harmonic number
                    
                    Ucompl = sys.real_to_compl(U);
                    
                    D1vec = [zeros(Hu+1,1) ; reshape((1:Hu).*ones(2*Hu+1,1),2*Hu*Hu+Hu,1)];
                    D2vec = [(0:1:Hu)' ; reshape((-Hu:1:Hu)'.*ones(1,Hu),2*Hu*Hu+Hu,1)];
                    
                    % Complex :
                    DUcompl = 1i*(omega1*D1vec + omega2*D2vec).*Ucompl;
                    % Real :
                    DU = sys.compl_to_real(DUcompl);
            end
        end
        
        function [DDU] = DD(sys,U,omega1sq,omega2sq,omega1omega2)
            %  Compute  D(D.U) only for QPHBM systems.
            Hu=floor(sqrt((size(U,1)-1)/4));   % harmonic number
            
            Ucompl = sys.real_to_compl(U);
            
            D1vec = [zeros(Hu+1,1) ; reshape((1:Hu).*ones(2*Hu+1,1),2*Hu*Hu+Hu,1)];
            D2vec = [(0:1:Hu)' ; reshape((-Hu:1:Hu)'.*ones(1,Hu),2*Hu*Hu+Hu,1)];
            
            % Complex :
            DDUcompl = -(omega1sq*(D1vec.*D1vec) + omega2sq*(D2vec.*D2vec) + (2*omega1omega2)*(D1vec.*D2vec)).*Ucompl;
            % Real :
            DDU = sys.compl_to_real(DDUcompl);
        end
            
        function R = Prod(sys,U,V,complex)
            % Compute the product of U and V using complex convolution
            switch sys.type
                case 'HBM'
                    U = full(U);
                    V = full(V);
                    
                    if nargin < 4
                        Ucompl = sys.real_to_compl(U,'full');
                        Vcompl = sys.real_to_compl(V,'full');
                    else
                        Ucompl = U;
                        Vcompl = V;
                    end
                    
                    Rcompl_tot = zeros(size(U));
                    for i=1:size(U,2)
                        Rcompl_tot(:,i) = conv2(Ucompl(:,i),Vcompl(:,i),'same');
                    end
                    Rcompl = Rcompl_tot(sys.H+1:end,:);
                    R = sys.compl_to_real(Rcompl);
                    
                case 'QPHBM'
                    Ucompl = sys.real_to_compl(U);
                    Vcompl = sys.real_to_compl(V);
                    
                    % Transform vectors into matrices of coefficients
                    matU = sys.get_mat_var(Ucompl);
                    matV = sys.get_mat_var(Vcompl);
                    
                    % Perform the central part of 2D convolution
                    matR = conv2(matU,matV,'same');
                    
                    % Create the vector version
                    Rcompl = matR(:);
                    Hu=floor(sqrt((size(U,1)-1)/4));   % harmonic number
                    Rcompl = Rcompl((2*Hu+1)*Hu+Hu+1:end);
                    R = sys.compl_to_real(Rcompl);
            end
        end
        
        function Bmat = B(sys,U)
            %  Compute  the matrix B from U (size 2H+1 )
            %
            switch sys.type
                case 'HBM'
                    %  toeplitz(c,r) generate a matrix from a colum c and a row r
                    %
                    Hu=(size(U,1)-1)/2   ;   % harmonic number
                    Uc= U(2:Hu+1)/2;           % Half of the cosines
                    Us= U(Hu+2:2*Hu+1)/2;       % Half of the sines
                    
                    LUc= [ U(1) ; Uc(1:end-1) ];    %   [ U0 ; Uc1/2 ; Uc2/2 ...]
                    LUs= [   0  ; Us(1:end-1) ];    %   [ 0  ; Us1/2 ; Us2/2 ...]
                    CC = toeplitz(LUc,LUc);
                    SS = toeplitz(-LUs,LUs);
                    
                    Cp = zeros(Hu,Hu) ;
                    Sp = zeros(Hu,Hu) ;
                    L0 = zeros(Hu,1) ;
                    VUc= [   0  ; Uc(end:-1:2) ];     %   [ 0  ; UcH/2 ; UcH-1/2 ...]
                    VUs= [   0  ; Us(end:-1:2) ];     %   [ 0  ; UsH/2 ; UsH-1/2 ...]
                    Cp(:,Hu:-1:1) = toeplitz(L0,VUc);
                    Sp(:,Hu:-1:1) = toeplitz(L0,VUs);
                    
                    Bmat(:,1)= U ;
                    
                    Bmat(1,2:Hu+1) = Uc ;
                    Bmat(1,Hu+2:2*Hu+1) = Us ;
                    
                    Bmat(2:Hu+1    ,2:Hu+1    )=CC + Cp;
                    Bmat(2:Hu+1    ,Hu+2:2*Hu+1)=SS + Sp;
                    Bmat(Hu+2:2*Hu+1,2:Hu+1    )=SS'+ Sp;
                    Bmat(Hu+2:2*Hu+1,Hu+2:2*Hu+1)=CC - Cp;
                case 'QPHBM'
                    Hu=floor(sqrt((size(U,1)-1)/4));   % harmonic number
                    
                    Ucompl = sys.real_to_compl(U);
                    Umat = sys.get_mat_var(Ucompl);
                    Bmtx = convmtx2(Umat,(2*Hu+1),(2*Hu+1));
                    
                    NT = Hu*(4*Hu+1);
                    Bmtx(end-NT+1:end,:)=[];
                    list_remove = [(1:NT)'; reshape(repmat(NT+(1:Hu)',1,(2*Hu+1))+(0:1:2*Hu)*(4*Hu+1),(2*Hu+1)*Hu,1) ; ...
                        reshape(repmat(NT+(1:Hu)'+(3*Hu+1),1,(2*Hu+1))+(0:1:2*Hu)*(4*Hu+1),(2*Hu+1)*Hu,1)];
                    Bmtx(list_remove,:)=[];
                    % Now Bmtx is such that : reshape(Bmtx*Vmat(:),2*H+1,2*H+1) = conv2(Umat,Vmat,'same');
                    nb_harmo = 2*Hu*(Hu+1);
                    B00 = Bmtx(nb_harmo+1,nb_harmo+1);
                    B0h = Bmtx(nb_harmo+1,1+nb_harmo+1:end);
                    Bh0 = Bmtx(1+nb_harmo+1:end,nb_harmo+1);
                    Bhmh= Bmtx(1+nb_harmo+1:end,nb_harmo:-1:1);
                    Bhh = Bmtx(1+nb_harmo+1:end,1+nb_harmo+1:end);
                    
                    % Featuring :
                    Bmat = [B00         , real(B0h)         , imag(B0h) ;
                        2*real(Bh0) , real(Bhh + Bhmh)  , imag(Bhh - Bhmh) ;
                        -2*imag(Bh0), -imag(Bhh + Bhmh)  , real(Bhh - Bhmh) ];
            end
        end
        
        function Ureal = compl_to_real(sys,Ucompl)
            % function Ureal = compl_to_real(sys,Ucompl)
            % Take a vector of complex harmonics  from 0 to H and give
            % the real corresponding Fourier series (constant, cos, sin)
            Ucos = 2*real(Ucompl(2:end,:));
            Usin = -2*imag(Ucompl(2:end,:));
            Ureal = real([Ucompl(1,:); Ucos ; Usin]);
        end
        
        function Ucompl = real_to_compl(sys,Ureal,shape)
            % function Ucompl = real_to_compl(sys,Ureal)
            % It is the inverse function of compl_to_real
            ind_midvec = (size(Ureal,1)+1)/2; % position of the last cosine.
            
            Uharm = 0.5*(Ureal(2:ind_midvec,:)-1i*Ureal(ind_midvec+1:end,:));
            if nargin < 3
                Ucompl = [Ureal(1,:);Uharm];
            else
                Ucompl = [flipud(conj(Uharm));Ureal(1,:);Uharm];
            end
        end
        
        function [ M ] = get_mat_var( sys, V  )
            %Transforms a vector representing indices :
            % (0,0);(0,1);..;(0,H);(1,-H);(1,H);(2,-H);......;(H,H)
            %into a (2H+1)*(2H+1)-matrix of indices between -H and H:
            % (-H,-H)   ...     (-H, H)
            %           ...
            % (H, -H)   ...     (H, H)
            Hv=floor(sqrt( (size(V,1)-1)/2));   % harmonic number
            
            if size(V,2) == 1
                V_0w2pos = V(1:Hv+1);
                v1 = V(Hv+2:end);
                
                Mpos = reshape(v1,[2*Hv+1 Hv]);
                
                M = zeros(2*Hv+1);
                
                M(:,1:Hv) = rot90(conj(Mpos),2);
                M(:,Hv+1) = [flipud(conj(V_0w2pos(2:end)));V_0w2pos];
                M(:,Hv+2:end) = Mpos;
            else
                V_0w2pos = V(1:Hv+1,:);
                v1 = V(Hv+2:end,:);
                
                Mpos = reshape(v1,[2*Hv+1 Hv size(V,2)]);
                
                M = zeros(2*Hv+1,2*Hv+1,size(V,2));
                
                M(:,1:Hv,:) = rot90(conj(Mpos),2);
                M(:,Hv+1,:) = [flipud(conj(V_0w2pos(2:end,:)));V_0w2pos];
                M(:,Hv+2:end,:) = Mpos;
            end
            
        end
        
    end
end


