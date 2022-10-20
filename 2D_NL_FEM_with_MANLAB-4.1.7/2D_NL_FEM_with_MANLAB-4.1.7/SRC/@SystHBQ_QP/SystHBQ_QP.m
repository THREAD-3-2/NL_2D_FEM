classdef SystHBQ_QP < Syst
    
    properties
        
        subtype = 'autonomous' ; % subtype of system: {['autonomous'],'forced','conservative?'} for later use.
        
        H;          % Number of harmonics
        nz_tot ;    % number of equation of the DAE
        nz ;        % number of main equations
        nz_aux=0;   % number of auxiliary equations
        
        pulsation='omega';  % for forced system, the continuation parameter is omega by default.
        % Otherwise, sys.pulsation is the value of the
        % pulsation of the forcing terms
        
        % Small operator on the Algebro-differential equations
        id =[]; jd =[]; vd =[];                % list for the sparse (order 2) tensor d
        id1=[]; jd1=[]; vd1=[];                % list for the sparse (order 2) tensor d1
        idd=[]; jdd=[]; vdd=[];                % list for the sparse (order 2) tensor dd
        ic0=[]; vc0=[];                        % list for the sparse (order 1) tensor c0
        ic1=[]; vc1=[];                        % list for the sparse (order 1) tensor c1
        ic2=[]; vc2=[];                        % list for the sparse (order 1) tensor c2
        il0=[]; jl0=[]; vl0=[];                % list for the sparse (order 2) tensor l0
        il1=[]; jl1=[]; vl1=[];                % list for the sparse (order 2) tensor l1
        iq =[]; jq =[]; kq =[]; vq=[];         % list for the sparse (order 3) tensor q
        
        iforce0=[]; hforce0=[]; vforce0=[];    % list for the sparse (order 2) tensor force
        iforce1=[]; hforce1=[]; vforce1=[];    % list for the sparse (order 2) tensor force
        iforce2=[]; hforce2=[]; vforce2=[];    % list for the sparse (order 2) tensor force
        
        idl=[]; jdl=[]; vdl=[];                % list for the sparse (order 2) tensor dl
        idq=[]; jdq=[]; kdq=[]; vdq=[];        % list for the sparse (order 3) tensor dq
                
        zi_phase1=1; % index of the state variable receiving the first phase condition
        zi_phase2=1; % index of the state variable receiving the second phase condition
    end
    
    methods
        
        function sys = SystHBQ_QP(nz,nz_aux,H,equations,point_display,global_display,parameters,subtype,writing)
            
            if nargin<8; subtype = 'autonomous'; end
            if nargin<9; writing = 'standard'; end
            
            coef_per_var = (H+1)*4*H+1;
            sys=sys@Syst('neq',nz*coef_per_var+2,'neq_aux',nz_aux*coef_per_var+5);
            
            sys.type = 'HBQ_QP';
            
            sys.parameters = parameters;
            sys.equations = equations;
            sys.point_display = point_display;
            sys.global_display = global_display;
            
            sys.nz = nz;
            sys.nz_aux = nz_aux;
            sys.nz_tot = nz + nz_aux;
            sys.H  = H ;    % harmonic number
            
            % Arclength only on the main variables
            sys.arclengthdef = sparse(1:sys.neq+1,ones(sys.neq+1,1),1,sys.ninc,1);
            
            % Creation of the little operators
            switch writing
                case 'vectorial'
                    sys = sys.get_operators_vec;
                case 'standard'
                    sys = sys.get_operators;
                otherwise
                    error('SystHBQ : Uncorrect writing type of the system. It should be either standard or vectorial.');
            end
            
            % subtype
            sys.subtype = subtype;
            
            % Pulsation of the oscillations for forced systems.
            sys.subtype = subtype;
            if strcmp(sys.subtype,'forced')
                try sys.pulsation = sys.parameters.pulsation;
                catch; sys.pulsation = 'omega'; disp('The focring pulsation is the continuation parameter.');
                end
            end
            
            % Residu
            sys.R = @R;
            
            %% Computation of the residue
            function [Rtot] = R(obj,Utot)
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
                        Rphase1= Ztot((coef_per_var-1)/2+1+(1+2*H),sys.zi_phase1);       % sin(omega1) de la variable zi_phase1
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
        %% All the function used to compute the residu and the jacobian from the little operators
        
        function [Ztot,omega1,omega2,lambda,varargout] = get_Ztot(obj,Utot)
            
            coef_per_var = (obj.H+1)*4*obj.H+1;
            Z      = reshape(Utot(1:obj.neq-2),coef_per_var,obj.nz);
            omega1 =Utot(obj.neq-1);
            omega2 =Utot(obj.neq);
            lambda=Utot(obj.neq+1) ;
            
            varargout{1} = Utot(end-4); % omega1sq =Utot(end-4);
            varargout{2} = Utot(end-3); % omega2sq =Utot(end-3);
            varargout{3} = Utot(end-2); % omega1omega2 =Utot(end-2);
            varargout{4} = Utot(end-1); % lambdaomega1 =Utot(end-1);
            varargout{5} = Utot(end);   % lambdaomega2 =Utot(end);
            Zaux = reshape(Utot(obj.neq+2:end-5),coef_per_var,obj.nz_aux);
            
            Ztot = [Z , Zaux];
        end
        
        function  [DU] = D(obj,U,omega1,omega2)
            %  Compute  D.U
            Hu=floor(sqrt((size(U,1)-1)/4));   % harmonic number
            
            Ucompl = obj.real_to_compl(U);
            
            D1vec = [zeros(Hu+1,1) ; reshape((1:Hu).*ones(2*Hu+1,1),2*Hu*Hu+Hu,1)];
            D2vec = [(0:1:Hu)' ; reshape((-Hu:1:Hu)'.*ones(1,Hu),2*Hu*Hu+Hu,1)];
            
            % Complex :
            DUcompl = 1i*(omega1*D1vec + omega2*D2vec).*Ucompl;
            % Real :
            DU = obj.compl_to_real(DUcompl);
        end
        
        function  [DDU] = DD(obj,U,omega1sq,omega2sq,omega1omega2)
            %  Compute  D(D.U)
            Hu=floor(sqrt((size(U,1)-1)/4));   % harmonic number
            
            Ucompl = obj.real_to_compl(U);
            
            D1vec = [zeros(Hu+1,1) ; reshape((1:Hu).*ones(2*Hu+1,1),2*Hu*Hu+Hu,1)];
            D2vec = [(0:1:Hu)' ; reshape((-Hu:1:Hu)'.*ones(1,Hu),2*Hu*Hu+Hu,1)];
            
            % Complex :
            DDUcompl = -(omega1sq*(D1vec.*D1vec) + omega2sq*(D2vec.*D2vec) + (2*omega1omega2)*(D1vec.*D2vec)).*Ucompl;
            % Real :
            DDU = obj.compl_to_real(DDUcompl);
        end
        
        function R = Prod(obj,U,V)
            % Compute the product of U and V using complex convolution
            
            U = full(U);
            V = full(V);
            
            Ucompl = obj.real_to_compl(U);
            Vcompl = obj.real_to_compl(V);
            
            % Transform vectors into matrices of coefficients
            matU = obj.get_mat_var(Ucompl);
            matV = obj.get_mat_var(Vcompl);
            
            % Perform the central part of 2D convolution
            matR = conv2(matU,matV,'same');
            
            % Create the vector version
            Rcompl = matR(:);
            Hu=floor(sqrt((size(U,1)-1)/4));   % harmonic number
            Rcompl = Rcompl((2*Hu+1)*Hu+Hu+1:end);
            R = obj.compl_to_real(Rcompl);
        end
        
        function Bmat = B(obj,U)
            % Compute  the matrix B from U (size 2H+1 )
            % With convmtx :
            Hu=floor(sqrt((size(U,1)-1)/4));   % harmonic number
            
            Ucompl = obj.real_to_compl(U);
            Umat = obj.get_mat_var(Ucompl);
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
        
        function Ureal = compl_to_real(obj,Ucompl)
            % function Ureal = compl_to_real(sys,Ucompl)
            % Take a vector of complex harmonics  from 0 to H and give
            % the real corresponding Fourier series (constant, cos, sin)
            
            Ucos = 2*real(Ucompl(2:end,:));
            Usin = -2*imag(Ucompl(2:end,:));
            Ureal = real([Ucompl(1,:); Ucos ; Usin]);
        end
        
        function Ucompl = real_to_compl(obj,Ureal)
            % function Ucompl = real_to_compl(sys,Ureal)
            % It is the inverse function of compl_to_real
            Hu=floor(sqrt( (size(Ureal,1)-1)/4));   % harmonic number
            Ucompl = [Ureal(1,:);0.5*Ureal(2:2*Hu*(Hu+1)+1,:)-0.5*1i*Ureal(2*Hu*(Hu+1)+2:end,:)];
        end
        
        function [ M ] = get_mat_var( obj, V  )
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


