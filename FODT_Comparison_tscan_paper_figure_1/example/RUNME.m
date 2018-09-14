path(path,'../MAT');
proc_surf_plots

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
