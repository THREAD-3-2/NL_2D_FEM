classdef SystHBQ < Syst
    
    properties
        
        H;  % Order of truncation of Fourier series
        
        subtype = 'autonomous' ; % subtype of system: {['autonomous'],'forced'}
        
        nz_tot ;    % number of equation of the DAE
        nz ;        % number of main equations
        nz_aux=0;   % number of auxiliary equations
        
        angfreq='omega';  % for forced system, the continuation parameter is omega by default.
        % Otherwise, sys.angfreq is the value of the
        % angular frequency of the forcing terms
        
        % Small operator on the Algebro-differential equations
        id =[]; jd =[]; vd= [];                % list for the sparse (order 2) tensor d
        id1=[]; jd1=[]; vd1=[];                % list for the sparse (order 2) tensor d1
        idd=[]; jdd=[]; vdd=[];                % list for the sparse (order 2) tensor dd
        ic0=[]; vc0=[];                        % list for the sparse (order 1) tensor c0
        ic1=[]; vc1=[];                        % list for the sparse (order 1) tensor c1
        ic2=[]; vc2=[];                        % list for the sparse (order 1) tensor c1
        il0=[]; jl0=[]; vl0=[];                % list for the sparse (order 2) tensor l0
        il1=[]; jl1=[]; vl1=[];                % list for the sparse (order 2) tensor l1
        iq =[]; jq =[]; kq =[]; vq=[];         % list for the sparse (order 3) tensor q
        uiq=[];i_iq=[];
        
        iforce0=[]; hforce0=[]; vforce0=[];       % list for the sparse (order 2) tensor force
        iforce1=[]; hforce1=[]; vforce1=[];       % list for the sparse (order 2) tensor force
        iforce2=[]; hforce2=[]; vforce2=[];       % list for the sparse (order 2) tensor force
        
        idl=[]; jdl=[]; vdl=[];                % list for the sparse (order 2) tensor dl
        idq=[]; jdq=[]; kdq=[]; vdq=[];        % list for the sparse (order 3) tensor dq
        
        zi_phase=1; % index of the state variable reciving the phase condition zi'(0)=0
        % or if the system is 'conservative', it receives zi(0) = lambda,
        % while other initial conditions are free.
    end
    
    methods
        
        function sys = SystHBQ(nz,nz_aux,H,equations,point_display,global_display,parameters,subtype,writing)
            
            sys=sys@Syst('neq',nz*(2*H+1)+1,'neq_aux',nz_aux*(2*H+1)+2);
            
            sys.type = 'HBQ';
            
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
            
            % subtype - bifurcation following
            if nargin > 7
                sys.subtype = subtype;
            end
            
            % Angular frequency of the oscillations for forced systems.
            if strcmp(sys.subtype,'forced')
                try sys.angfreq = sys.parameters.angfreq;
                catch; sys.angfreq = 'omega';
                end
            end
            
            % Creation of the little operators
            switch writing
                case 'vectorial'
                    sys = sys.get_operators_vec;
                case 'standard'
                    sys = sys.get_operators;
                otherwise
                    error('SystHBQ : Uncorrect writing type of the system. It should be either standard or vectorial.');
            end
            
            % Specific treatment of the products to avoid a long loop in Fpnl.m
            [sys.uiq,sys.i_iq,~] = unique(sys.iq); % Select one product per equation to allow a vectorial allocation
            
            
            
            % Residu
            sys.R = @R;
            
            function [Rtot] = R(sys,Uf)
                %  Compute  R(U)
                H    = sys.H;   % Harmonic number
                nz_tot   = sys.nz_tot;  % number of equation of the DAE system
                DHp1 = 2*H+1 ;
                
                % Extract informations from Utot
                [Zf,omega,lambda,omega2,lambdaomega] = sys.get_Ztot(Uf);
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
                    Rfourier(:,sys.il0(i))= Rfourier(:,sys.il0(i)) + sys.vl0(i)*Zf(:,sys.jl0(i)) ;
                end
                
                for i=1:numel(sys.il1)
                    Rfourier(:,sys.il1(i))= Rfourier(:,sys.il1(i)) + (sys.vl1(i)*lambda)*Zf(:,sys.jl1(i)) ;
                end
                
                for i=1:numel(sys.id)
                    Rfourier(:,sys.id(i))= Rfourier(:,sys.id(i)) + (sys.vd(i)*omega)*sys.D(Zf(:,sys.jd(i))) ;
                end
                
                for i=1:numel(sys.id1)
                    Rfourier(:,sys.id1(i))= Rfourier(:,sys.id1(i)) + (sys.vd1(i)*lambdaomega)*sys.D(Zf(:,sys.jd1(i))) ;
                end
                
                for i=1:numel(sys.idd)
                    Rfourier(:,sys.idd(i))= Rfourier(:,sys.idd(i)) + (sys.vdd(i)*omega2)*sys.D(sys.D(Zf(:,sys.jdd(i)))) ;
                end
                
                for i=1:numel(sys.iq)
                    Rfourier(:,sys.iq(i))= Rfourier(:,sys.iq(i)) + sys.vq(i)*( sys.Prod(Zf(:,sys.jq(i)),Zf(:,sys.kq(i))) ) ;
                end
                
                Ztot_t0 = sum(Zf(1:H+1,:),1);
                
                % phase equation
                switch sys.subtype
                    case 'autonomous'
                        Rphase= (1:H)*Zf(H+2:2*H+1,sys.zi_phase);
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
                        Rfourier(:,sys.idl(i))= Rfourier(:,sys.idl(i)) + sys.vdl(i)*sys.D(Zf(:,sys.jdl(i))) ;
                    end
                    
                    for i=1:numel(sys.idq)
                        Rfourier(:,sys.idq(i))= Rfourier(:,sys.idq(i)) + sys.vdq(i)*( sys.Prod(Zf(:,sys.jdq(i)),sys.D(Zf(:,sys.kdq(i)))) ) ;
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
            
            
            
        end
        %% All the function used to compute the residu and the jacobian from the little operators
        function [Ztot,omega,lambda,omega2,lambdaomega] = get_Ztot(sys,Utot)
            
            DHp1 = sys.H*2+1;
            Z      = reshape(Utot(1:sys.neq-1),DHp1,sys.nz);
            omega  = Utot(sys.neq);
            lambda = Utot(sys.neq+1) ;
            
            omega2 = Utot(end-1);
            lambdaomega =Utot(end);
            Zaux   = reshape(Utot(sys.neq+2:end-2),DHp1,sys.nz_aux);
            
            Ztot = [Z Zaux];
            
        end
        
        function  DU = D(obj,U,complex)
            %  Compute  D.U
            Hu=(size(U,1)-1)/2   ;   % harmonic number
            if nargin < 3
                DU = [ zeros(1,size(U,2)) ; (1:Hu)'.*U(Hu+2:end,:) ; (-(1:Hu))'.*U(2:Hu+1,:) ];
            else
                DU = 1i*((-Hu:Hu)'.*U);
            end
        end
        
        function R = Prod(obj,U,V,complex)
            % Compute the product of U and V using complex convolution
            
            U = full(U);
            V = full(V);
            
            if nargin < 4
                Ucompl = obj.real_to_compl(U,'full');
                Vcompl = obj.real_to_compl(V,'full');
            else
                Ucompl = U;
                Vcompl = V;
            end
            
            Rcompl_tot = zeros(size(U));
            for i=1:size(U,2)
                Rcompl_tot(:,i) = conv2(Ucompl(:,i),Vcompl(:,i),'same');
            end
            Rcompl = Rcompl_tot(obj.H+1:end,:);
            R = obj.compl_to_real(Rcompl);
        end
        
        function Bmat = B(obj,U)
            %  Compute  the matrix B from U (size 2H+1 )
            %
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
        end
        
        function Ureal = compl_to_real(obj,Ucompl)
            % function Ureal = compl_to_real(sys,Ucompl)
            % Take a vector of complex harmonics  from 0 to H and give
            % the real corresponding Fourier series (constant, cos, sin)
            Ucos = 2*real(Ucompl(2:end,:));
            Usin = -2*imag(Ucompl(2:end,:));
            Ureal = real([Ucompl(1,:); Ucos ; Usin]);
        end
        
        function Ucompl = real_to_compl(obj,Ureal,shape)
            % function Ucompl = real_to_compl(sys,Ureal)
            % It is the inverse function of compl_to_real
            Hu=(size(Ureal,1)-1)/2   ;   % harmonic number
            Uharm = 0.5*(Ureal(2:Hu+1,:)-1i*Ureal(Hu+2:end,:));
            if nargin < 3
                Ucompl = [Ureal(1,:);Uharm];
            else
                Ucompl = [flipud(conj(Uharm));Ureal(1,:);Uharm];
            end
        end
        
        
    end
end


