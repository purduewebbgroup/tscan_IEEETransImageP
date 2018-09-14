function out=mysurf(X,Y,f, colorz, map_flag)
figure;
set(gcf, 'DefaultAxesFontSize',30);
surf(X,Y,f');
view(0,90);
shading flat;
xlabel('x (cm)');
ylabel('y (cm)');
caxis(colorz);
m=[255:-1:0]'/255;
if map_flag==1
  colormap([m m m]);
else
  colormap('gray');
end 
pos1=get(gca,'Position');
set(gca,'Position', pos1+[.017 .038 -.07 .00])

cb=colorbar;
pos2=get(cb, 'Position');
set(cb,'Position', pos2+[-.02 .00 -0.00 0])
