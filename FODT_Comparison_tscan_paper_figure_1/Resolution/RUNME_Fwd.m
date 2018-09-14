path(path,'../MAT');
%proc_surf_plots

clc; clear all; close all

FS = 48;

Nq = 6;
Nm = 6;
N = Nq*Nm;
foo= 'meas';
Data = loadmeas(foo,N);

Real_Data = Data(:,1);
Imag_Data = Data(:,2);

Mag_Data = (Real_Data.^2 + Imag_Data.^2).^0.5;
Phase_Data = atan(Imag_Data./Real_Data);

figure
plot(Real_Data);
xlabel('Data Index','fontweight','bold','fontsize',FS);
ylabel('Real(Y)','fontweight','bold','fontsize',FS);
title('Fwd Solution Data','fontweight','bold','fontsize',FS);
set(gca,'fontweight','bold','fontsize',FS);
axis tight

figure
plot(Imag_Data);
xlabel('Data Index','fontweight','bold','fontsize',FS);
ylabel('Imag(Y)','fontweight','bold','fontsize',FS);
title('Fwd Solution Data','fontweight','bold','fontsize',FS);
set(gca,'fontweight','bold','fontsize',FS);
axis tight

figure
plot(Mag_Data);
xlabel('Data Index','fontweight','bold','fontsize',FS);
ylabel('Mag(Y)','fontweight','bold','fontsize',FS);
title('Fwd Solution Data','fontweight','bold','fontsize',FS);
set(gca,'fontweight','bold','fontsize',FS);
axis tight

figure
plot(Phase_Data);
xlabel('Data Index','fontweight','bold','fontsize',FS);
ylabel('Phase(Y)','fontweight','bold','fontsize',FS);
title('Fwd Solution Data','fontweight','bold','fontsize',FS);
set(gca,'fontweight','bold','fontsize',FS);
axis tight

save Y_nn_d1 Data