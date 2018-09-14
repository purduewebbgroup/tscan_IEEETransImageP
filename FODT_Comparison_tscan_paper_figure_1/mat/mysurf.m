function out=mysurf4(X,Y,f, colorz, map1, map2, xlbl, ylbl, tlbl)
% function out=mysurf4(X,Y,f, colorz, map1, map2, xlbl, ylbl)

figure;
set(gcf, 'DefaultAxesFontSize',30);
surf(X,Y,f');
Xrange=max(X)-min(X);
Yrange=max(Y)-min(Y);
set(gca,'XLim', [min(X)-Xrange/45 max(X)+Xrange/45]);
set(gca,'YLim', [min(Y)-Yrange/45 max(Y)+Yrange/45]);
daspect([1 1 1]);
view(0,90);
shading flat;
xlabel(xlbl);
ylabel(ylbl);
title(tlbl);
caxis(colorz);
for k=0:255;
  m(k+1,1:3)=map1(1,1:3)+(map2(1,1:3)-map1(1,1:3))*k/255;
end
colormap(m);
%pos1=get(gca,'Position');
%set(gca,'Position', pos1+[.017 .038 -.07 .00])
set(gca,'Box','on')
set(gca,'LineWidth',[1.0])

cb=colorbar;
pos2=get(cb, 'Position');
set(cb,'Position', pos2+[-.02 .00 -0.00 0])


