%% Brian Bentz 9/2015
% Generate phys file

% Sources on bottom and Detectors on top

clc; clear all; close all;
FS = 48;

%% OUTPUT
fn = sprintf('phys');
fid = fopen(fn,'w+');

%% Scater (Background)
scale = 1;

% Medium (Brain) 
n = 1.33;
mus = 20/scale;           % cm^-1
mua = 0.2/scale;        % cm^-1
zo = 1/mus;
l_mfp = 1/mus;

D = 1/(3*(mus + mua))  


%% Physical Parameters

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

delta_x = (xmax-xmin)/(Ni-1) % These should equal l_mpf
delta_y = (ymax-ymin)/(Nj-1)
delta_z = (zmax-zmin)/(Nl-1)

ls = 32*l_mfp;               % medium height along z in mm
ld = 32*l_mfp;               % medium width along x in mm
ly = 32*l_mfp;               % medium width along y in mm


%% Define Src/Det Positions
% for speed set K and M to number of src and det
K = 6;
M = 6;

% 4 Source 4 Detectors
position_det = [8*l_mfp 16*l_mfp 0; ...
                  0 16*l_mfp 24*l_mfp; ...
                  24*l_mfp 16*l_mfp ls; ...
                  ld 16*l_mfp 8*l_mfp;
                  24*l_mfp 0 16*l_mfp;
                  8*l_mfp ly 16*l_mfp];
position_src = [24*l_mfp 16*l_mfp 0; ...
                  0 16*l_mfp 8*l_mfp; ...
                  8*l_mfp 16*l_mfp ls; ...
                  ld 16*l_mfp 24*l_mfp;
                  8*l_mfp 0 16*l_mfp;
                  24*l_mfp ly 16*l_mfp];

% Sources and Detectors on each side
figure
for index = 1:6   
    x = position_src(index,1);
    y = position_src(index,2);
    z = position_src(index,3);
    hold on;plot3(x,y,z,'o','linewidth',1,'MarkerEdgeColor','w','MarkerFaceColor','r');

    x = position_det(index,1);
    y = position_det(index,2);
    z = position_det(index,3);
    hold on;plot3(x,y,z,'o','linewidth',1,'MarkerEdgeColor','w','MarkerFaceColor','b');
end
st = sprintf('FODT Sources (red) and Detectors (Blue)');
title(st,'interpreter','none','fontweight','bold','fontsize',FS);
xlabel('Width: x [cm]','fontweight','bold','fontsize',FS);
ylabel('Length: y [cm]','fontweight','bold','fontsize',FS);
zlabel('Height: z [cm]','fontweight','bold','fontsize',FS);
set(gca,'fontweight','bold','fontsize',FS);
axis([xmin xmax ymin ymax zmin zmax])
daspect([1 1 1]);

Nq = K;
Nm = M;


%% Output File

K = Nq;  % number of sources
M = Nm;  % number of detectors

c = 2.998*10^10; % spead of light in cm/s
v = c/n;  % velocity: c/1.3333

% output to file
fprintf(fid,'# Physical Parameters\n');
fprintf(fid,'physical {\n');
fprintf(fid,' xmin  %g\n',xmin);
fprintf(fid,' xmax   %g\n',xmax);
fprintf(fid,' ymin  %g\n',ymin);
fprintf(fid,' ymax   %g\n',ymax);
fprintf(fid,' zmin  %g\n',zmin);
fprintf(fid,' zmax   %g\n',zmax);
fprintf(fid,'\n');
fprintf(fid,'  Ni  %g\n',Ni);
fprintf(fid,'  Nj  %g\n',Nj);
fprintf(fid,'  Nl  %g\n',Nl);
fprintf(fid,'\n');
fprintf(fid,'  K  %g\n',K);
fprintf(fid,'  M  %g\n',M);
fprintf(fid,'\n');
fprintf(fid,'  v   %g\n',v); % delete + sign
fprintf(fid,'}\n');
fprintf(fid,'\n');
fprintf(fid,'\n');

%% Sources
omega = 4.927988476219283e+08;
beta = 1.0;
calir = 1.0;
calii = 0.0;

% output to file
fprintf(fid,'#Sources\n');
fprintf(fid,'\n');

for q = 1:Nq
    fprintf(fid,'source {\n');
    fprintf(fid,'#%g\n',q-1);
    fprintf(fid,'  position %g %g %g\n',position_src(q,1),position_src(q,2),position_src(q,3));
    fprintf(fid,'  omega  %.15e\n',omega);
    fprintf(fid,'  beta  %.1f\n',beta);
    fprintf(fid,'  calir  %.1f\n',calir);
    fprintf(fid,'  calii  %.1f\n',calii);
    fprintf(fid,'}\n');
end

fprintf(fid,'\n');
fprintf(fid,'\n');
fprintf(fid,'\n');


%% Detectors
omega = 4.927988476219283e+08;
beta = 1.0;
calir = 1.0;
calii = 0.0;

% output to file
fprintf(fid,'#detectors\n');

for m = 1:Nm
    fprintf(fid,'detector {\n');
    fprintf(fid,'#%g\n',m-1);
    fprintf(fid,'  position %g %g %g\n',position_det(m,1),position_det(m,2),position_det(m,3));
    fprintf(fid,'  omega  %.15e\n',omega);
    fprintf(fid,'  calir  %.1f\n',calir);
    fprintf(fid,'  calii  %.1f\n',calii);
    fprintf(fid,'}\n');
end

fprintf(fid,'\n');
fprintf(fid,'\n');


%% Connections
conect = [0:1:Nm];

fprintf(fid,'#Connections\n');
fprintf(fid,'connections {\n');

for q = 1:Nq
    fprintf(fid,'%d [',q-1);
    for m = 1:Nm
        fprintf(fid,' %d',m-1);
    end
    fprintf(fid,' ]\n');    
end

fprintf(fid,'}\n');


