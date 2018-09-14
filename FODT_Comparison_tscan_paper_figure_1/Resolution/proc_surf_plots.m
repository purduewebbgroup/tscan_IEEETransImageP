get_img;                                                                               

t=6.7;
T=max(size(t));

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

xmin = -l_mfp;
xmax = (Ni-2)*l_mfp;
ymin = -l_mfp;
ymax = (Nj-2)*l_mfp;
zmin = -l_mfp;
zmax = (Nl-2)*l_mfp;

delta_x = (xmax-xmin)/(Ni-1); % These should equal l_mpf
delta_y = (ymax-ymin)/(Nj-1);
delta_z = (zmax-zmin)/(Nl-1);
                                                                               
% X=-4.18:8.36/32:4.18;
% Y=-4.18:8.36/32:4.18;
% Z=-0.18:(5.88+.18)/16:5.88;

X=xmin:delta_x:xmax;
Y=ymin:delta_y:ymax;
Z=zmin:delta_z:zmax;

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
