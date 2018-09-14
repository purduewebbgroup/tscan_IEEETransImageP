% X=-4.18:8.36/32:4.18;
% Y=-4.18:8.36/32:4.18;
% Z=-0.18:(5.88+.18)/16:5.88;

z1_top=zeros(33,33);
z2_top=zeros(33,33);
z1_bott=zeros(33,33);
z2_bott=zeros(33,33);

% z1_top=zeros(65,65);
% z2_top=zeros(65,65);
% z1_bott=zeros(65,65);
% z2_bott=zeros(65,65);

foo=['OUTPUT/fluohat'];
x=loadmu1(foo,'fluohat',33,33,33);
%x=loadmu1(foo,'fluohat',65,65,33);
[x1,trash]=shiftdim(x(1,:,:,:));
[x2,trash]=shiftdim(x(2,:,:,:));

%   z1_top(:,:)=x1(:,:,12);
%   z2_top(:,:)=x2(:,:,12);
%   z1_bott(:,:)=x1(:,:,7);
%   z2_bott(:,:)=x2(:,:,7);

%   z1_top(:,:)=x1(:,:,9);
%   z2_top(:,:)=x2(:,:,9);
%   z1_bott(:,:)=x1(:,:,9);
%   z2_bott(:,:)=x2(:,:,9);

z1_top(:,:)=x1(:,:,16);
z2_top(:,:)=x2(:,:,16);
z1_bott(:,:)=x1(:,:,16);
z2_bott(:,:)=x2(:,:,16);

