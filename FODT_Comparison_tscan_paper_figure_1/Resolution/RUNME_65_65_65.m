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
Ni = 65; %(2^6 + 1)
Nj = 65;
Nl = 65;

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

% Units of l_mfp
X = X./l_mfp;
Y = Y./l_mfp;
Z = Z./l_mfp;


%% Make Actual Positions 
y1_actual=zeros(Ni,Nj,Nl);

x = 0.746875;
y = 0.800000;
z = 0.746875;

dx = round(x/delta_x);
dy = round(y/delta_y);
dz = round(z/delta_z);

y1_actual(dx,dy,dz) = 0.1;  % mm^-1

x = 0.746875;
y = 0.800000;
z = 0.853125;

dx = round(x/delta_x);
dy = round(y/delta_y);
dz = round(z/delta_z);

y1_actual(dx,dy,dz) = 0.15;  % mm^-1

x = 0.853125;
y = 0.800000;
z = 0.746875;

dx = round(x/delta_x);
dy = round(y/delta_y);
dz = round(z/delta_z);

y1_actual(dx,dy,dz) = 0.075;  % mm^-1

x = 0.853125;
y = 0.800000;
z = 0.853125;

dx = round(x/delta_x);
dy = round(y/delta_y);
dz = round(z/delta_z);

y1_actual(dx,dy,dz) = 0.05;  % mm^-1


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

ya_actual(:,:)=y1_actual(:,30,:);

y1_bott2(:,:)=x1(:,30,:);
y2_bott2(:,:)=x2(:,30,:);
y1_bott(:,:)=x1(:,31,:);
y2_bott(:,:)=x2(:,31,:);
y1_top(:,:)=x1(:,32,:);
y2_top(:,:)=x2(:,32,:);


%% Plot surface
F=max(max(max(x2)));
F1=max(max(max(x2)));

Ta = max(max(ya_actual))*1.2;
T1 = max(max(y2_top));
T2 = max(max(y2_bott));
T3 = max(max(y2_bott2));

mysurf2(X,Y,ya_actual,[0,0.15],[1 1 1], [0 0 0], 'x (l*)', 'z (l*)','\eta');
title('$\eta$','Interpreter','LaTex')

mysurf(X,Y,y2_top,[0 T1],[1 1 1], [0 0 0], 'x (l*)', 'z (l*)','\eta');
mysurf(X,Y,y2_bott,[0 T2],[1 1 1], [0 0 0], 'x (l*)', 'z (l*)','\eta');
mysurf(X,Y,y2_bott2,[0 T3],[1 1 1], [0 0 0], 'x (l*)', 'z (l*)','\eta');


% Plot Isoimage
figure
iso_value = 0.03; % cm^-1
st = sprintf('Eta = %.2f',iso_value);
niceplot2(X,Y,Z,x2,iso_value,[1 0 1],30,st)
title('$\hat{\eta} = 0.003$','Interpreter','LaTex')
view(0,0)

figure
iso_value = 0.003; % mm^-1
st = sprintf('Eta = %.2f',iso_value);
niceplot2(X,Y,Z,y1_actual,iso_value,[1 0 1],30,st)
title('$\hat{\eta} = 0.003$','Interpreter','LaTex')
view(0,0)


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
