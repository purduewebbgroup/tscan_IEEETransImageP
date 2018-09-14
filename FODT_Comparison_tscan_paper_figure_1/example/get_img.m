X=-4.18:8.36/32:4.18;
Y=-4.18:8.36/32:4.18;
Z=-0.18:(5.88+.18)/16:5.88;

z1_top=zeros(33,33);
z2_top=zeros(33,33);
z1_bott=zeros(33,33);
z2_bott=zeros(33,33);

  foo=['OUTPUT/fluohat'];
  x=loadmu1(foo,'fluohat',33,33,17);
  [x1,trash]=shiftdim(x(1,:,:,:));
  [x2,trash]=shiftdim(x(2,:,:,:));
  z1_top(:,:)=x1(:,:,12);
  z2_top(:,:)=x2(:,:,12);
  z1_bott(:,:)=x1(:,:,7);
  z2_bott(:,:)=x2(:,:,7);

