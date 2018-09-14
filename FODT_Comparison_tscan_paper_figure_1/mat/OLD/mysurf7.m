function out=mysurf7(X,Y,f, colorz, map1, map2)
%function out=mysurf7(X,Y,f, colorz, map1, map2)
figure;
%set(gcf, 'DefaultAxesFontSize',25);
set(gcf, 'DefaultAxesFontSize',30);
surf(X,Y,f');
view(0,90);
shading flat;
xlabel('x (cm)');
ylabel('y (cm)');
caxis(colorz);
for k=0:255;
%  m(k+1,1:3)=map1(1,1:3)+((map2(1,1:3)-map1(1,1:3))*k/255);
  m(k+1,1:3)=map1(1,1:3)+(map2(1,1:3)-map1(1,1:3))*((k/255).^2);
end
m
%m=m.^.75;
%m=m.^.45;
colormap(m);
set(gca,'XLimMode','manual');
set(gca,'XLim',[min(X)-.4 max(X)+.4]);
set(gca,'YLimMode','manual');
set(gca,'YLim',[min(Y)-.4 max(Y)+.4]);
pos1=get(gca,'Position');
set(gca,'Position', pos1+[.017 .038 -.07 .00])
set(gca,'Box','on');
set(gca,'LineWidth',1.75);
daspect([1 1 .6])

%cb=colorbar;
%pos2=get(cb, 'Position');
%set(cb,'Position', pos2+[-.02 .00 -0.00 0])
%set(cb,'LineWidth',1.75);
