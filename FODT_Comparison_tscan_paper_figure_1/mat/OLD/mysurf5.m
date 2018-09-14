function out=mysurf4(X,Y,f, colorz, map1, map2)
figure;
% set(gcf, 'DefaultAxesFontSize',20);
 set(gcf, 'DefaultAxesFontSize',30);
surf(X,Y,f');
view(0,90);
shading flat;
xlabel('x (cm)');
ylabel('y (cm)');
caxis(colorz);
for k=0:255;
  m(k+1,1:3)=map1(1,1:3)+(map2(1,1:3)-map1(1,1:3))*k/255;
end
colormap(m);
pos1=get(gca,'Position');
set(gca,'Position', pos1+[.017 .038 -.07 .00])

cb=colorbar;
pos2=get(cb, 'Position');
set(cb,'Position', pos2+[-.02 .00 -0.00 0])
