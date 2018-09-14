path(path,'../MAT');

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


%% Make Actual Positions 
y1_actual=zeros(Ni,Nj,Nl);

x = 0.746875;
y = 0.800000;
z = 0.746875;

dx = round(x/delta_x);
dy = round(y/delta_y);
dz = round(z/delta_z);

y1_actual(dx,dy,dz) = 1.0;

x = 0.746875;
y = 0.800000;
z = 0.853125;

dx = round(x/delta_x);
dy = round(y/delta_y);
dz = round(z/delta_z);

y1_actual(dx,dy,dz) = 1.5;

x = 0.853125;
y = 0.800000;
z = 0.746875;

dx = round(x/delta_x);
dy = round(y/delta_y);
dz = round(z/delta_z);

y1_actual(dx,dy,dz) = 0.75;

x = 0.853125;
y = 0.800000;
z = 0.853125;

dx = round(x/delta_x);
dy = round(y/delta_y);
dz = round(z/delta_z);

y1_actual(dx,dy,dz) = 0.5;


%% Get Data
y1_top=zeros(Ni,Nl);
y2_top=zeros(Ni,Nl);
y1_bott=zeros(Ni,Nl);
y2_bott=zeros(Ni,Nl);

ya_actual = zeros(Ni,Nl);

foo=['OUTPUT/fluohat'];
x=loadmu1(foo,'fluohat',Ni,Nj,Nl);
[x1,trash]=shiftdim(x(1,:,:,:));
[x2,trash]=shiftdim(x(2,:,:,:));

ya_actual(:,:)=y1_actual(:,15,:);

y1_bott2(:,:)=x1(:,15,:);
y2_bott2(:,:)=x2(:,15,:);
y1_bott(:,:)=x1(:,16,:);
y2_bott(:,:)=x2(:,16,:);
y1_top(:,:)=x1(:,17,:);
y2_top(:,:)=x2(:,17,:);


%% Plot surface
F=max(max(max(x2)));
F1=max(max(max(x2)));

Ta = max(max(ya_actual));
T1 = max(max(y2_top))/3;
T2 = max(max(y2_bott))/3;
T3 = max(max(y2_bott2))/3;

mysurf(X,Y,ya_actual,[0,Ta],[1 1 1], [0 0 0], 'x (cm)', 'z (cm)','Eta 1');

mysurf(X,Y,y2_top,[0 T1],[1 1 1], [0 0 0], 'x (cm)', 'z (cm)','Eta 1');
mysurf(X,Y,y2_bott,[0 T2],[1 1 1], [0 0 0], 'x (cm)', 'z (cm)','Eta 1');
mysurf(X,Y,y2_bott2,[0 T3],[1 1 1], [0 0 0], 'x (cm)', 'z (cm)','Eta 1');


% Plot Isoimage
figure
iso_value = max(max(max(x2)))/3;
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
