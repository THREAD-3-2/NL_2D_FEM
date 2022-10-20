% function Hand = plotfloquetmultHBM(sys,Section,bifpara_str)
%    plots the Floquet multipliers at the end of the Section
%    as a function of the bifurcation parameter
%    specified by bifpara_str ( = 'omega' or 'lambda').
%
%    It also displays the Floquet at the bifurcation points.
%
%    By L. Guillot and O. Thomas / Nov. 2018

function Hand = plotfloquetmultHBM(sys,Sec_Diag,bifpara_str,Comparison_method)

if nargin < 4
    Comparison_method = 'imag';
end

subplot(2,2,1)
hold on
subplot(2,2,3)
hold on
subplot(2,2,[2 4])
hold on
ang = linspace(0,2*pi,500);
plot(cos(ang),sin(ang),'k-')

% Color for successive plots.
Col = dispcolorML(sys.nz);

switch bifpara_str
    case 'omega'
        ind_bifpara = sys.neq;
    case 'lambda'
        ind_bifpara = sys.neq+1;
end

if isa(Sec_Diag,'CheckPoint')
    Hand = plotsectionHBM(sys,Sec_Diag);
else
    cellfun(@(sec)plotsectionHBM(sys,sec),Sec_Diag(1:end-1),'uniformoutput',0);
    Hand = plotsectionHBM(sys,Sec_Diag{end});
end


subplot(2,2,1)
xlabel(['\' bifpara_str])
ylabel('Modulus of the Floquet multipliers');
subplot(2,2,3)
xlabel(['\' bifpara_str])
ylabel('Argument of the Floquet multipliers');
subplot(2,2,[2 4])
xlabel('Real part of the Floquet multipliers');
ylabel('Imaginary part of the Floquet multipliers');


    function Hand = plotsectionHBM(sys,Section)
        
        Hand = [];
        
        floq_exp = Section.Eigen_end.values;
        switch Comparison_method
            case 'imag'
                [~,ind] = sort(imag(floq_exp));
            case 'real'
                [~,ind] = sort(real(floq_exp));
        end
        floq_exp = floq_exp(ind);
        
        Omega_end = Section.Uend(sys.neq);
        floq_mult = exp(2*pi/Omega_end*floq_exp);
        
        %%% Test if there is a bifurcation in the Section
        if ~strcmp(Section.Eigen.type,'nothing')
            floq_exp_stab = Section.Eigen.values;
            switch Comparison_method
                case 'imag'
                    [~,ind] = sort(imag(floq_exp_stab));
                case 'real'
                    [~,ind] = sort(real(floq_exp_stab));
            end
            floq_exp_stab = floq_exp_stab(ind);
            
            Omega_end = Section.Uend(sys.neq);
            floq_mult_stab = exp(2*pi/Omega_end*floq_exp_stab);
        end
        
        Nb = length(floq_exp);
        
        for i=1:Nb
            cdisp = Col(i,:);
            
            % Real part
            subplot(2,2,1)
            hand = plot(Section.Uend(ind_bifpara),abs(floq_mult(i)),'.','Color',cdisp);
            Hand = [Hand;hand];
            
            % Imaginary part
            subplot(2,2,3)
            hand = plot(Section.Uend(ind_bifpara),angle(floq_mult(i)),'.','Color',cdisp);
            Hand = [Hand;hand];
            
            % Plot the multipliers in the complex plane.
            subplot(2,2,[2 4])
            hand = plot(real(floq_mult(i)),imag(floq_mult(i)),'.','Color',cdisp);
            Hand = [Hand;hand];
            
            %%% Test if there is a bifurcation in the Section
            if ~strcmp(Section.Eigen.type,'nothing')
                % Real part
                subplot(2,2,1)
                hand = plot(Section.Ustab(ind_bifpara),abs(floq_mult_stab(i)),'p','Color',cdisp);
                Hand = [Hand;hand];
                
                % Imaginary part
                subplot(2,2,3)
                hand = plot(Section.Ustab(ind_bifpara),angle(floq_mult_stab(i)),'p','Color',cdisp);
                Hand = [Hand;hand];
                
                % Plot the multipliers in the complex plane.
                subplot(2,2,[2 4])
                hand = plot(real(floq_mult_stab(i)),imag(floq_mult_stab(i)),'p','Color',cdisp);
                Hand = [Hand;hand];
            end
        end
    end

end