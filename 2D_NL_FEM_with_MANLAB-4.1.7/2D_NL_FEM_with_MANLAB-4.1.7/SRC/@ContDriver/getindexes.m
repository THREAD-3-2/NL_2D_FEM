function [whatfound,idxsection,Apoint] = getindexes(obj)
figure(2);
[x,y]=ginput(1);  % get the clicked point
d = inf;
nbp = length(obj.ChckPoint);

% for each checkpoint, look for the distance between all the point used
% in the ProjBif Diag (regular and bif) and the cliked point

for i=1:nbp
    X = get(obj.ChckPoint{i},'X');   % X(nbpts,ncurve)
    Y = get(obj.ChckPoint{i},'Y');   % Y(nbpts,ncurve)
    D = sqrt(abs(X - x).^2 + abs(Y - y).^2);     % D(nbpts,ncurve) distance of each point
    % of the projBifDiag with the clicked point
    [dist,idx] =  min(D) ;  %  dist(ncurve) minimal value for each column of D
    %  idx(ncurve) index in each column
    [mindist,k]=min(dist);  % mindist : min of D , found on row k of D
    if (mindist < d)
        d = mindist;
        whatfound  = 'regular';
        idxsection = i ;
        idxpoint   = idx(k);
        idxcurve   = k ;
        A          = get(obj.ChckPoint{idxsection},'A'); % A(nbpts)
        Apoint     = A(idxpoint);
    end
    
    chkp_status = get(obj.ChckPoint{i},'BifStatus');
    
    if strcmp(chkp_status{1},'nothing')==0
        Xbif = get(obj.ChckPoint{i},'Xbif');   % Xbif(ncurve)
        Ybif = get(obj.ChckPoint{i},'Ybif');   % Ybif(ncurve)
        Xstab = get(obj.ChckPoint{i},'Xstab'); % Xstab(ncurve)
        Ystab = get(obj.ChckPoint{i},'Ystab'); % Ystab(ncurve)
        Dbif = sqrt(abs(Xbif - x).^2 + abs(Ybif - y).^2);  % Dbif(ncurve) distance of each point
        % bif point of the projBifDiag with the clicked point
        Dstab = sqrt(abs(Xstab - x).^2 + abs(Ystab - y).^2);
        % Dbif is divided by two to increase the probability to choose a SB point.
        [mindbif,kbif]=min(Dbif/2);  % mindist : min of Dbif , found on row k of Dbif
        [mindstab,kstab]=min(Dstab);  % mindist : min of Dstab , found on row k of Dstab
        if mindbif < min([mindstab,d])
            d = mindbif;
            whatfound  = 'simplebif';
            idxsection = i ;
            idxpoint   = 0;
            idxcurve   = kbif ;
            Apoint     = 0 ;
        elseif mindstab < min([mindbif , d])
            d = mindstab;
            whatfound = chkp_status{1};
            idxsection = i ;
            idxpoint   = 0;
            idxcurve   = kstab ;
            Apoint     = get(obj.ChckPoint{i},'Astab');
        end
    end
end




switch whatfound
    case 'regular'
        disp([' Regular point selected : Section ', num2str(idxsection), ...
            ' Point ' , num2str(idxpoint),' Curve ' , num2str(idxcurve) ]);
        
    case 'simplebif'
        disp(['Branching point selected : Section ', num2str(idxsection), ...
            ' Curve ' , num2str(idxcurve)]) ;
        
    case 'Hopf'
        disp(['Hopf point selected : Section ', num2str(idxsection), ...
            ' Curve ' , num2str(idxcurve)]) ;
        
    case 'PD'
        disp(['Period-doubling point selected : Section ', num2str(idxsection), ...
            ' Curve ' , num2str(idxcurve)]) ;
        
    case 'SB'
        disp(['Simple Bifurcation point selected : Section ', num2str(idxsection), ...
            ' Curve ' , num2str(idxcurve)]) ;
        
    case 'NS'
        disp(['Neimarck-Sacker point selected : Section ', num2str(idxsection), ...
            ' Curve ' , num2str(idxcurve)]) ;

end

end

