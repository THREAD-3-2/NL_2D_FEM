function [Diag] = calcdiagUpp(sys,Diagram)
% function [Diag] = calcdiagUpp(sys,Diagram)
%
% Puts the point representation of all the sections in on big matrix
% DiagUpp. If there are only one output argument, it returns the
% informations extracted from the Diagram in the structure "Diag".
%
% Stabinfo is a line vector of length size(DiagUpp,2) which contains
% stability informations : 0 if unstable solution, 1 if stable solution,
% double('B') = 66 if simple bifurcation point (branching points or folds), 
% double('H') = 72 if Hopf bifurcation point, 
% double('P') = 80 if period doubling , 
% double('N') = 78 if Neimark-Sacker bifurcation,
% NaN if the bifurcation could not be determined correctly.
%
% If a jump is detected between two solution branches, a bifurcation point
% with associated character ' ' is plotted. double(' ') = 32.
%
% The inverse function of "double" is "char" in this context : char(78)='N'.

if ~iscell(Diagram)
    Diagg{1} = Diagram; % Transform a single Section in a Diagram with a single Section.
    Diagram=Diagg;
end

nb_Section = length(Diagram);
nb_pts = size(Diagram{1}.Upp,2);
nb_pts_tot = zeros(length(Diagram),1);

DiagUpp = zeros(sys.ninc,nb_Section*nb_pts);
Stabinfo = zeros(1,nb_Section*nb_pts);

UUend=Diagram{1}.Upp(:,1);
%disp(['end0 ',num2str(UUend(2))])
nbjump=0;

for kk = 1:nb_Section
    Sectionkk = Diagram{kk};

    % Stability test 
    if Sectionkk.Ustab == 0 
        nb_pts_tot(kk) = nb_pts;
        % Discrete representation of the solution-branch and stability.
        Uppkk = Sectionkk.Upp;
        Stabinfo(sum(nb_pts_tot(1:kk-1))+(1:nb_pts_tot(kk))+nbjump) = strcmp(Sectionkk.StabStatus{1},'stable');
    else
        nb_pts_tot(kk) = nb_pts+1;
        ind_change = Sectionkk.ind_change;
        % Discrete representation of the solution-branch.
        Uppkk = [Sectionkk.Upp(:,1:ind_change-1), Sectionkk.Ustab, Sectionkk.Upp(:,ind_change:end)];
        
        % Stability information at the beginning of the branch
        Stabinfo( sum(nb_pts_tot(1:kk-1))+(1:ind_change-1) +nbjump) = strcmp(Sectionkk.StabStatus{1},'stable');
        % Type of bifurcation : 'B', 'H', 'PD' or 'NS' 
        % [simple bifurcation, Hopf, period doubling, Neimark-Sacker]
        Stabinfo( sum(nb_pts_tot(1:kk-1)) + ind_change +nbjump) = double(Sectionkk.Eigen.type(1));
        % Stability information at the end of the branch
        Stabinfo( sum(nb_pts_tot(1:kk-1)) + (ind_change+1:nb_pts_tot(kk)) +nbjump ) = strcmp(Sectionkk.StabStatus{2},'stable');
    end
    II=sum(nb_pts_tot(1:kk-1))+(1:nb_pts_tot(kk))+nbjump; % vector of index of the present section
    UU0=Uppkk(:,1);
    testjump= norm(UUend-UU0)>0; % test for a jump in the diagram
    DiagUpp(:,II)= Uppkk;
    if testjump
        nbjump=nbjump+1;
        disp(['Jump detected before section ',num2str(kk)])
        DiagUpp=[DiagUpp(:,1:(II(1)-1)) NaN*ones(sys.ninc,1) DiagUpp(:,II(1):end)];
        Stabinfo=[Stabinfo(1:(II(1)-1)) , double(' ') , Stabinfo(II(1):end)];
    end
    UUend=Uppkk(:,end);
end

% Replace all the value double('n') = 110 by NaN, meaning that the
% detection of the bifurcation failed.
Stabinfo(Stabinfo == 110) = NaN;

% Array of index(ices) of the stability change(s)
change = find(Stabinfo > 1);

if nargout == 1
    Diag.DiagUpp = DiagUpp;
    Diag.Stabinfo = Stabinfo;
    Diag.change = change;
else
    Diag = DiagUpp;
end

