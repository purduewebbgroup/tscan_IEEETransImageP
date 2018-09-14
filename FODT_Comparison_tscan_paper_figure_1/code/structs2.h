typedef struct fcn_arglist_type{
  prior_param *priorparam;
  mtype ****x;
  int u;
  int i;
  int j;
  int l;
  int R;
  mtype **theta1;
  mtype **theta2; 
  mtype *neighbor;
  mtype *freqs;
} arglist; 

