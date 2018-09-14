%get_img;  

close all;

path(path,'../MAT');

z1_top=zeros(33,33);
z2_top=zeros(33,33);
z1_bott=zeros(33,33);
z2_bott=zeros(33,33);

% z1_top=zeros(65,65);
% z2_top=zeros(65,65);
% z1_bott=zeros(65,65);
% z2_bott=zeros(65,65);

%foo=['OUTPUT/fluohat'];
%x=loadmu1(foo,'fluohat',33,33,17);
%x=loadmu1(foo,'fluohat',65,65,33);

foo=['fluo'];
x=loadmu1(foo,'fluo',33,33,17);

[x1,trash]=shiftdim(x(1,:,:,:));
[x2,trash]=shiftdim(x(2,:,:,:));
                                                                               
% X=-4.18:8.36/32:4.18;
% Y=-4.18:8.36/32:4.18;
% Z=-0.18:(5.88+.18)/16:5.88;

X=-0.518:1.036/32:0.518;
Y=-0.618:1.236/32:0.618;
Z=-0.05:(0.45+.05)/16:0.45;

% X=-0.518:1.036/64:0.518;
% Y=-0.618:1.236/64:0.618;
% Z=-0.05:(0.45+.05)/32:0.45;

F=max(max(max(x2)));
F1=max(max(max(x2)))/3;
close all;
for dex=9:9    
    z1_top(:,:)=x1(:,:,dex);
    z2_top(:,:)=x2(:,:,dex);
    z1_bott(:,:)=x1(:,:,dex);
    z2_bott(:,:)=x2(:,:,dex);    
    mysurf(X,Y,z2_bott,[0 F],[1 1 1], [0 0 0], 'x (cm)', 'y (cm)','Eta 1');
    print('-deps', ['bottom_blob.eps'])
    %close;
    mysurf(X,Y,z2_top,[0 F1],[1 1 1], [0 0 0], 'x (cm)', 'y (cm)','Eta 1');
    print('-deps', ['top_blob.eps'])
    %close;
end

% Plot Isoimage
figure
iso_value = max(max(max(x2)))/3;
st = sprintf('Eta = %d',iso_value);
niceplot(X,Y,Z,x2,iso_value,[1 0 1],30,st)
