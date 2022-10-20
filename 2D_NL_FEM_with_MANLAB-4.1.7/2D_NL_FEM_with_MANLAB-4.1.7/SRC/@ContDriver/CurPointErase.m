function h =CurPointErase(obj,params,h)
% erase the current point 
figure(2); hold on;
s=size(params.dispvars);
for i=1:s(1)
   delete( h{i} );
end

h=cell(0);
