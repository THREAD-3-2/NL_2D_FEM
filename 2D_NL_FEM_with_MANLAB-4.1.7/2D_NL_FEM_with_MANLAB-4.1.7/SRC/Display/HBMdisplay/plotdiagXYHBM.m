% function Hand = plotdiagXYHBM(sys,Diag,X,Y)
%    plots each line of the variables X,Y in a 2D plot.
%    X and Y are matrix of size Nb*nbpts_diag where Nb is the number of
%    curve wanted and nbpts_diag is the number of points along the diagram
%    (obtained with the function calcdiagUpp).
%
%    It also displays the bifurcation points:
%    - B for a simple bifurcation (saddle-node, pitchfork...)
%    - PD for period doubling bifurcation
%    - NS for Neimark-Sacker bifurcation
%
%
%    By L. Guillot and O. Thomas / Nov. 2018 - Mars 2019


function Hand = plotdiagXYHBM(sys,Diag,X,Y)

if ~isstruct(Diag)
    [Diag] = calcdiagUpp(sys,Diag);
end

DiagUpp = Diag.DiagUpp;
Stabinfo = Diag.Stabinfo;
change = Diag.change;

NbX=size(X,1);
NbY=size(Y,1);
if NbX ~= NbY
    error('plotdiagXYHBM.m : X ans Y must be of the same size.');
end
Nb = NbX;

Hand=[];

% Color for successive plots.
Col = dispcolorML2(Nb);


strstab = '-';
strunstab = ':';

for ii=1:Nb
    
    cdisp = Col(ii,:);
    Xdisp = X(ii,:);
    Ydisp = Y(ii,:);
    if numel(change) == 0
        
        if Stabinfo(1) == 1; strplot = strstab; else; strplot = strunstab; end
        hand = plot(Xdisp,Ydisp,strplot,'Color',cdisp,'linewidth',1);
        Hand = [Hand;hand];
        
    else
        
        if numel(change) == 0
            if Stabinfo(1) == 1; strplot = strstab; else; strplot = strunstab; end
            
            hand = plot(Xdisp,Ydisp,strplot,'Color',cdisp,'linewidth',1);
            Hand = [Hand;hand];
            
        else
            for istab=1:numel(change)
                ind_change = change(istab);
                if istab == 1; indplot = 1:ind_change; else; indplot = change(istab-1):ind_change; end
                if Stabinfo(indplot(2)) == 1; strplot = strstab; else; strplot = strunstab; end
                
                hand = plot(Xdisp(indplot),Ydisp(indplot),strplot,'Color',cdisp,'linewidth',1);
                Hand = [Hand;hand];
                hand = plot(Xdisp(ind_change),Ydisp(ind_change),'.','Color',cdisp);
                Hand = [Hand;hand];
                
%text(bifparastab,Ystab,char(Stabinfo(ind_change)),'margin',3,'verticalalignment','Bottom')
                box on;
            end
            
            if Stabinfo(ind_change+1) == 1; strplot = strstab; else; strplot = strunstab; end
            hand = plot(Xdisp(ind_change:end),Ydisp(ind_change:end),strplot,'Color',cdisp,'linewidth',1);
            Hand = [Hand;hand];
            
        end
    end
end




end