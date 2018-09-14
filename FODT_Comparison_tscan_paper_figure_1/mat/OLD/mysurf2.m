function mysurf2(XX,YY,x1);

surf(XX,YY,x1'); 
daspect([1 1 .4]); 
view(0,90) 
m=(255:-1:0)/255;
mymap=[m' m' m']
colormap(mymap);
