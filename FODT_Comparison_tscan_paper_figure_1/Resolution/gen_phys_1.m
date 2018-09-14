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
l_mfp = 1/mus;

D = 1/(3*(mus + mua));  % alph = 1: 0.0167 (cm)
                       % alph = 2: 0.0333 (cm)

%% Physical Parameters

Ni = 33; %(2^5 + 1)
Nj = 33;
Nl = 33;

% Ni = 65; % (2^6 + 1)
% Nj = 65;
% Nl = 65;

Bounded = 1;
xmin = -l_mfp*Bounded;
xmax = 32*l_mfp + l_mfp*Bounded;
ymin = -l_mfp*Bounded;
ymax = 32*l_mfp + l_mfp*Bounded;
zmin = -l_mfp*Bounded;
zmax = 32*l_mfp + l_mfp*Bounded;

delta_x = (xmax-xmin)/(Ni-1); % These should equal l_mpf
delta_y = (ymax-ymin)/(Nj-1);
delta_z = (zmax-zmin)/(Nl-1);


%% Define Src/Det Positions
% for speed set K and M to number of src and det
K = 118;
M = 120;
position_src = ones(K,3);  % xpos,ypos,zpos
position_det = ones(M,3);  % xpos,ypos,zpos

Nx = 5; % number of src and det along x dimension
Ny = 5;
Nz = 5;

% del_sx_xy_xz should all be > zo
del_sx = (xmax-xmin)/(Nx-1)
del_sy = (ymax-ymin)/(Ny-1)
del_sz = (zmax-zmin)/(Nz-1)
zo

% Sources and Detectors on each side
figure
Nq = 0;
Nm = 0;
% src/det on top/bot
count = 0;
for x = (xmin + del_sx : del_sx : xmax - del_sx)
    for y = (ymin + del_sy : del_sy : ymax - del_sy)
        count = count + 1;        
        if mod(count,2) == 0 % alternate src det
            z = zmin + zo; % bottom
            Nq = Nq + 1; % number of sources        
            position_src(Nq,1) = x;
            position_src(Nq,2) = y;
            position_src(Nq,3) = z;
            hold on;plot3(x,y,z,'o','linewidth',1,'MarkerEdgeColor','w','MarkerFaceColor','r');
            
            z = zmax - zo; % top
            Nq = Nq + 1; % number of sources        
            position_src(Nq,1) = x;
            position_src(Nq,2) = y;
            position_src(Nq,3) = z;
            hold on;plot3(x,y,z,'o','linewidth',1,'MarkerEdgeColor','w','MarkerFaceColor','r');
        else            
            z = zmin + zo; % bottom
            Nm = Nm + 1; % number of detectors
            position_det(Nm,1) = x;
            position_det(Nm,2) = y;
            position_det(Nm,3) = z; 
            hold on;plot3(x,y,z,'o','linewidth',1,'MarkerEdgeColor','w','MarkerFaceColor','b');
            
            z = zmax - zo; % top
            Nm = Nm + 1; % number of detectors
            position_det(Nm,1) = x;
            position_det(Nm,2) = y;
            position_det(Nm,3) = z;
            hold on;plot3(x,y,z,'o','linewidth',1,'MarkerEdgeColor','w','MarkerFaceColor','b');
        end
    end
end


% src/det on xz planes (2 sides)
% count = 0; % remove?
for x = (xmin + del_sx : del_sx : xmax - del_sx)
    for z = (zmin + del_sz : del_sz : zmax - del_sz)
        count = count + 1;        
        if mod(count,2) == 0 % alternate src det
            y = ymin + zo; % ymin side
            Nq = Nq + 1; % number of sources        
            position_src(Nq,1) = x;
            position_src(Nq,2) = y;
            position_src(Nq,3) = z;
            hold on;plot3(x,y,z,'o','linewidth',1,'MarkerEdgeColor','w','MarkerFaceColor','r');
            
            y = ymax - zo; % ymax side
            Nq = Nq + 1; % number of sources        
            position_src(Nq,1) = x;
            position_src(Nq,2) = y;
            position_src(Nq,3) = z;
            hold on;plot3(x,y,z,'o','linewidth',1,'MarkerEdgeColor','w','MarkerFaceColor','r');
        else            
            y = ymin + zo; % ymin side
            Nm = Nm + 1; % number of detectors
            position_det(Nm,1) = x;
            position_det(Nm,2) = y;
            position_det(Nm,3) = z; 
            hold on;plot3(x,y,z,'o','linewidth',1,'MarkerEdgeColor','w','MarkerFaceColor','b');
            
            y = ymax - zo; % ymax side
            Nm = Nm + 1; % number of detectors
            position_det(Nm,1) = x;
            position_det(Nm,2) = y;
            position_det(Nm,3) = z;
            hold on;plot3(x,y,z,'o','linewidth',1,'MarkerEdgeColor','w','MarkerFaceColor','b');
        end
    end
end


% src/det on yz planes (2 sides)
%count = 0; % remove?
for y = (ymin + del_sy : del_sy : ymax - del_sy)
    for z = (zmin + del_sz : del_sz : zmax - del_sz)
        count = count + 1;        
        if mod(count,2) == 0 % alternate src det
            x = xmin + zo; % xmin side
            Nq = Nq + 1; % number of sources        
            position_src(Nq,1) = x;
            position_src(Nq,2) = y;
            position_src(Nq,3) = z;
            hold on;plot3(x,y,z,'o','linewidth',1,'MarkerEdgeColor','w','MarkerFaceColor','r');
            
            x = xmax - zo; % ymax side
            Nq = Nq + 1; % number of sources        
            position_src(Nq,1) = x;
            position_src(Nq,2) = y;
            position_src(Nq,3) = z;
            hold on;plot3(x,y,z,'o','linewidth',1,'MarkerEdgeColor','w','MarkerFaceColor','r');
        else            
            x = xmin + zo; % xmin side
            Nm = Nm + 1; % number of detectors
            position_det(Nm,1) = x;
            position_det(Nm,2) = y;
            position_det(Nm,3) = z; 
            hold on;plot3(x,y,z,'o','linewidth',1,'MarkerEdgeColor','w','MarkerFaceColor','b');
            
            x = xmax - zo; % ymax side
            Nm = Nm + 1; % number of detectors
            position_det(Nm,1) = x;
            position_det(Nm,2) = y;
            position_det(Nm,3) = z;
            hold on;plot3(x,y,z,'o','linewidth',1,'MarkerEdgeColor','w','MarkerFaceColor','b');
        end
    end
end
st = sprintf('FODT Sources (red) and Detectors (Blue) with alpha = %d',alpha);
title(st,'interpreter','none','fontweight','bold','fontsize',FS);
xlabel('Width: x [cm]','fontweight','bold','fontsize',FS);
ylabel('Length: y [cm]','fontweight','bold','fontsize',FS);
zlabel('Height: z [cm]','fontweight','bold','fontsize',FS);
set(gca,'fontweight','bold','fontsize',FS);
axis([xmin xmax ymin ymax zmin zmax])
daspect([1 1 1]);

Nq
Nm


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


