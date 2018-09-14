X=-4.18:8.36/32:4.18;
Y=-4.18:8.36/32:4.18;
Z=-0.18:(5.88+.18)/16:5.88;

z1_top=zeros(33,33,21);
z2_top=zeros(33,33,21);
z1_bott=zeros(33,33,21);
z2_bott=zeros(33,33,21);

for t=1:21
  t
  foo=['fluohat' num2str(t-1) ];
  x=loadmu1(foo,'fluohat',33,33,17);
  [x1,trash]=shiftdim(x(1,:,:,:));
  [x2,trash]=shiftdim(x(2,:,:,:));
  z1_top(:,:,t)=x1(:,:,12);
  z2_top(:,:,t)=x2(:,:,12);
  z1_bott(:,:,t)=x1(:,:,7);
  z2_bott(:,:,t)=x2(:,:,7);
end

