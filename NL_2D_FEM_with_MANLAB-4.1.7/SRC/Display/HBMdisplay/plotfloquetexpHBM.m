% function Hand = plotfloquetexpHBM(sys,Section,bifpara_str)
%    plots the Floquet exponents at the end of the Section
%    as a function of the bifurcation parameter
%    specified by bifpara_str ( = 'omega' or 'lambda').
%
%    It also displays the Floquet at the bifurcation points.
%
%    By L. Guillot and O. Thomas / Nov. 2018 - Mars 2019

function Hand = plotfloquetexpHBM(sys,Sec_Diag,bifpara_str,Comparison_method)

if nargin < 4
    Comparison_method = 'imag';
end

subplot(2,2,1)
hold on
subplot(2,2,3)
hold on
subplot(2,2,[2 4])
hold on

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
ylabel('Real part of the Floquet exponents');
subplot(2,2,3)
xlabel(['\' bifpara_str])
ylabel('Imaginary part of the Floquet exponents');
subplot(2,2,[2 4])
xlabel('Real part of the Floquet exponents');
ylabel('Imaginary part of the Floquet exponents');


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
        end
        
        Nb = length(floq_exp);
        
        for i=1:Nb
            cdisp = Col(i,:);
            
            % Real part
            subplot(2,2,1)
            hand = plot(Section.Uend(ind_bifpara),real(floq_exp(i)),'.','Color',cdisp);
            Hand = [Hand;hand];
            
            % Imaginary part
            subplot(2,2,3)
            hand = plot(Section.Uend(ind_bifpara),imag(floq_exp(i)),'.','Color',cdisp);
            Hand = [Hand;hand];
            
            % Plot the exponents in the complex plane.
            subplot(2,2,[2 4])
            hand = plot(real(floq_exp(i)),imag(floq_exp(i)),'.','Color',cdisp);
            Hand = [Hand;hand];
            
            %%% Test if there is a bifurcation in the Section
            if ~strcmp(Section.Eigen.type,'nothing')
                % Real part
                subplot(2,2,1)
                hand = plot(Section.Ustab(ind_bifpara),real(floq_exp_stab(i)),'p','Color',cdisp);
                Hand = [Hand;hand];
                
                % Imaginary part
                subplot(2,2,3)
                hand = plot(Section.Ustab(ind_bifpara),imag(floq_exp_stab(i)),'p','Color',cdisp);
                Hand = [Hand;hand];
                
                % Plot the exponents in the complex plane.
                subplot(2,2,[2 4])
                hand = plot(real(floq_exp_stab(i)),imag(floq_exp_stab(i)),'p','Color',cdisp);
                Hand = [Hand;hand];
            end
        end
    end

end