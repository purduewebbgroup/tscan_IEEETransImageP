%path(path,'../MAT');
path(path,'../mat');

close all;

%% Scater (Background)
n = 1.333;
alpha = 1;
mus = 20/alpha;     % cm^-1
mua = 0.02/alpha; % cm^-1
zo = 1/mus;
l_mfp = 1/mus;

D = 1/(3*(mus + mua));  % alph = 1: 0.0167 (cm)
                       % alph = 2: 0.0333 (cm)

%% Geometry
Ni = 33; %(2^5 + 1)
Nj = 33;
Nl = 33;

Bounded = 1;
xmin = -l_mfp*Bounded;
xmax = 32*l_mfp + l_mfp*Bounded;
ymin = -l_mfp*Bounded;
ymax = 32*l_mfp + l_mfp*Bounded;
zmin = -l_mfp*Bounded;
zmax = 32*l_mfp + l_mfp*Bounded;

delta_x = (xmax-xmin)/(Ni-1);
delta_y = (ymax-ymin)/(Nj-1);
delta_z = (zmax-zmin)/(Nl-1);
                                                                              
X=xmin:delta_x:xmax;
Y=ymin:delta_y:ymax;
Z=zmin:delta_z:zmax;


%% Get Data
y1_top=zeros(33,33);
y2_top=zeros(33,33);
y1_bott=zeros(33,33);
y2_bott=zeros(33,33);

% z1_top=zeros(65,65);
% z2_top=zeros(65,65);
% z1_bott=zeros(65,65);
% z2_bott=zeros(65,65);

foo=['OUTPUT/fluohat'];
x=loadmu1(foo,'fluohat',33,33,33);
%x=loadmu1(foo,'fluohat',65,65,33);
[x1,trash]=shiftdim(x(1,:,:,:));
[x2,trash]=shiftdim(x(2,:,:,:));

y1_top(:,:)=x1(:,16,:);
y2_top(:,:)=x2(:,16,:);
y1_bott(:,:)=x1(:,16,:);
y2_bott(:,:)=x2(:,16,:);


%% Plot surface
F=max(max(max(x2)));
F1=max(max(max(x2)))/3;
close all;
for dex=1:T
  mysurf(X,Y,z2_bott(:,:,1+(dex-1)*4),[0 F],[1 1 1], [0 0 0], 'x (cm)', 'y (cm)','Eta 1');
  %print('-deps', ['bottom_blob.eps'])
  %close;
  mysurf(X,Y,z2_top(:,:,1+(dex-1)),[0 F1],[1 1 1], [0 0 0], 'x (cm)', 'y (cm)','Eta 1');
  %print('-deps', ['top_blob.eps'])
  %close;
end

% Plot Isoimage
figure
iso_value = max(max(max(x2)))/10;
st = sprintf('Eta = %d',iso_value);
niceplot(X,Y,Z,x2,iso_value,[1 0 1],30,st)


%% Plot Alpha and Cost
load -ascii OUTPUT/RESULT
iter=RESULT(:,1);
cost=RESULT(:,2);
alpha=RESULT(:,3);
figure
plot(iter,cost);
xlabel('Iteration number');
ylabel('Cost function');
title('Cost function vs iterations');
figure
plot(iter,alpha);
xlabel('Iteration number');
ylabel('\alpha estimate');
title('\alpha estimate vs iterations');
