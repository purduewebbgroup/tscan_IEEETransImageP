function niceplot(x,y,z,v,c,colorz,fsize,tlbl)
%function niceplot(x,y,z,v,c,colorz,fsize);

[sizex sizey sizez]=size(v)

v2=zeros(sizey, sizex, sizez);
for i=1:sizez
  v2(:,:,i)=v(:,:,i).';
end


p = patch(isosurface(x, y, z, v2, c));
isonormals(x,y,z,v2,p)
%set(p, 'FaceColor', [.8 .8 .8] , 'EdgeColor', 'none');
%set(p, 'FaceColor', [0 .65 1] , 'EdgeColor', 'none');
%set(p, 'FaceColor', [0 .55 1] , 'EdgeColor', 'none');
set(p, 'FaceColor', colorz , 'EdgeColor', 'none');
daspect([1 1 1])

% set(gca, 'FontSize', 20);
 set(gca, 'FontSize', fsize);
% set(gca, 'FontSize', 30);
% set(gca, 'FontWeight', 'bold');
%  set(gca, 'LineWidth', 1.5);
  set(gca, 'LineWidth', 1.75);
% set(gca, 'LineWidth', 2);

xmin=min(min(min(x)));
ymin=min(min(min(y)));
zmin=min(min(min(z)));

xmax=max(max(max(x)));
ymax=max(max(max(y)));
zmax=max(max(max(z)));


axis([xmin xmax ymin ymax zmin zmax]);

xlabel('x (l*)');
ylabel('y (l*)');
zlabel('z (l*)');
title(tlbl)

view(3)
%view(52.5,30)

camlight 
lighting phong
set(gca, 'Box', 'On');
