function Col = dispcolorML(nb)

% Permet d'avoir un jeu de 5 couleur pas trop mal pour les traces
%
% ex: plot(x,y,'Color',col(ii,:))
    %Col=[0 0.7 0;0 0.4 0;0 0.5 0.5;0 0.2 1;0.6 0 0.6;0.9 0 0.9;0.7 0.2 0;1 0 0;1 0.4 0;1 0.7 0];

if nargin==1
    Col = dispcolorML2(nb);
%     Col = [.6 .9 .9 ; ...
%            .7  1 .5 ; ...
%            1  1   0 ;.5 1   0 ;.7  1  0 ; ...
%            .5 .5  1 ;.7 .5  1 ; ...
%            1  .5 .5 ;.5 .5 .5 ;.7 .5 .5 ; ...
%            1  .5  0 ;.5 .5  0 ;.7 .5  0 ; ...
%            1   0  1 ;.5  0  1 ;.7  0  1 ; ...
%            1   0 .5 ;.5  0 .5 ;.7  0 .5 ];
%     Col = Col(1+mod(0:5:(7*nb-1),19),:);
else
    Col = dispcolorML2(8);
    %Col=[0 0.25 0.65; 0.3 0.5 0;0.7 0 0.3;0.9 0.4 0;0.4 0.2 0.5];
end

Col=repmat(Col,100,1);

end