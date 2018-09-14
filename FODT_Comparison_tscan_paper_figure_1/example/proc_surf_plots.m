get_img;                                                                               


X=-4.18:8.36/32:4.18;
Y=-4.18:8.36/32:4.18;
Z=-0.18:(5.88+.18)/16:5.88;

t=6.7;
T=max(size(t));

close all;
                                                                               
X=-4.18:8.36/32:4.18;
Y=-4.18:8.36/32:4.18;
Z=-0.18:(5.88+.18)/16:5.88;

F=1.2870;
F1=0.9601;
  close all;
for dex=1:T
  mysurf(X,Y,z2_bott(:,:,1+(dex-1)*4),[0 F],[1 1 1], [0 0 0], 'x (cm)', 'y (cm)');
  print('-deps', ['bottom_blob.eps'])
  %close;
  mysurf(X,Y,z2_top(:,:,1+(dex-1)),[0 F1],[1 1 1], [0 0 0], 'x (cm)', 'y (cm)');
  print('-deps', ['top_blob.eps'])
  %close;
end
