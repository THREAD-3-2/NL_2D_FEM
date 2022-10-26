function h =CurPointDisp(obj,params,sys)
% plot the current point (with a square marks) and the traveling direction on the proj bif diagram

figure(2); hold on;
s = axis; lxe = abs(s(2)-s(1));lye = abs(s(4)-s(3)); %

h=cell(0);

dispvars=params.dispvars; s=size(dispvars);
U0now=obj.CurPoint.U0now; Utnow=obj.CurPoint.Utnow;


for i=1:s(1)
    idxX = dispvars(i,1); idxY = dispvars(i,2);
    u0x = U0now(idxX); u0y = U0now(idxY);
    utx = Utnow(idxX); uty = Utnow(idxY);
    
    
    lt = sqrt((utx/lxe)^2+(uty/lye)^2)*10;
    utx=utx/lt; uty=uty/lt; vtx = -uty*lxe/lye; vty =  utx*lye/lxe;
    
    h{i}=plot(u0x,u0y,'sk',[u0x,u0x+utx],[u0y, u0y+uty],'-k', ...
        [ u0x+0.9*utx+0.1*(vtx), u0x+utx, u0x+0.9*utx-0.1*(vtx)], ...
        [ u0y+0.9*uty+0.1*(vty), u0y+uty, u0y+0.9*uty-0.1*(vty)],'k');
end
end
