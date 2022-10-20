function Hand = plotphasediagHBM(sys,U,Idisp,Colordisp)
%    Hand = plotphasediagHBM(sys,U,Idisp,Colordisp) plots a phase diagram
%    of the periodic solution computed by the Manlab/HBM method.
%
%    U is the vector of unknowns of the final algebraic system solved
%    by Manlab.
%
%    Idisp contains the indices of the entries of the initial physical
%    unknown vector u(t) to be displayed. It must be of length 2 or 3.
%
%    It returns in Hand a column vector of handles to lineseries objects,
%    one handle per plotted line.
%
%    Example: Hand=plotphasediagHBM(sys,U,[1 2],1)
%
%    By L. Guillot and O. Thomas / Nov. 2018


Hand = [];

if nargin == 3
    Colordisp = 'b';%[.5 0 .9];%'r'
end

switch sys.type
    case {'HBQ','HBM'}
        H = sys.H;
        Z=sys.get_Ztot(U);
        
        time=linspace(0,1,100*H);
        time = time(:);
        nt=length(time);
        
        VectCos=cos(2*pi*time*(1:H));
        VectSin=sin(2*pi*time*(1:H));
        
        MAT=[ ones(nt,1) , VectCos , VectSin ];  % MAT(nt,2*H+1)
        Utime=MAT*Z(:,Idisp);  % Ut(nt,nz)
        
        if numel(Idisp) == 1
            Utime(:,2) = MAT*sys.D(Z(:,Idisp));
            hand=plot(Utime(:,1),Utime(:,2),'Color',Colordisp);
            Hand=[Hand;hand];
            xlabel(['z_{' num2str(Idisp) '}']);
            ylabel(['z_{' num2str(Idisp) "}'"]);
        
        elseif numel(Idisp) == 2
            hand=plot(Utime(:,1),Utime(:,2),'Color',Colordisp);
            Hand=[Hand;hand];
            xlabel(['z_{' num2str(Idisp(1)) '}']);
            ylabel(['z_{' num2str(Idisp(2)) '}']);
        elseif numel(Idisp) == 3
            hand=plot3(Utime(:,1),Utime(:,2),Utime(:,3),'Color',Colordisp);
            Hand=[Hand;hand];
            xlabel(['z_{' num2str(Idisp(1)) '}']);
            ylabel(['z_{' num2str(Idisp(2)) '}']);
            zlabel(['z_{' num2str(Idisp(3)) '}']);
            view(3);
        end
        
    case {'HBQ_QP','QPHBM'}
        disp('calcphasediagHBM : Quasi-periodic solutions not available yet.');
end


end

