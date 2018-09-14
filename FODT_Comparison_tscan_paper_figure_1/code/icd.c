#include <stdio.h>
#include <math.h>
#include "defs.h"

mtype obj_fn_fluo(
  phys_param *physparam,
  prior_param *prior_tau,
  prior_param *prior_gamma,
  prior_param *prior_D,
  prior_param *prior_mua,
  mtype ****fluo_lin,
  mtype ****mu_m,
  mtype alpha
)
{
  int u,i,j,l,n,S;
  mtype x_curr;
  mtype *neighbor;
  mtype prior[4];
  prior_param *priorparam;

  for(u=0; u<2; u++){
    if(u==0){
      priorparam=prior_mua;
    }
    else{
      priorparam=prior_D;
    }

    neighbor = (mtype *) malloc(sizeof(mtype)*priorparam->Nneighbor);

    prior[u]=0.0;
    for(i=1; i<physparam->Ni-1; i++)
    for(j=1; j<physparam->Nj-1; j++)
    for(l=1; l<physparam->Nl-1; l++){
      x_curr=mu_m[u][i][j][l];
      get_neighbor(u,i,j,l, priorparam, mu_m, neighbor);
      for (n=0; n<priorparam->Nneighbor; n++)
        prior[u] += priorparam->b[n]
               * pow( fabs(x_curr-neighbor[n]), priorparam->p ) ;
    }
    prior[u]= prior[u]/pow(priorparam->sigma,priorparam->p)/priorparam->p;
    free(neighbor);
  }

  for(u=0; u<2; u++){
    if(u==0){
      priorparam=prior_tau;
    }
    else{
      priorparam=prior_gamma;
    }

    neighbor = (mtype *) malloc(sizeof(mtype)*priorparam->Nneighbor);

    prior[u+2]=0.0;
    for(i=1; i<physparam->Ni-1; i++)
    for(j=1; j<physparam->Nj-1; j++)
    for(l=1; l<physparam->Nl-1; l++){
      x_curr=fluo_lin[u][i][j][l];
      get_neighbor(u,i,j,l, priorparam, fluo_lin, neighbor);
      for (n=0; n<priorparam->Nneighbor; n++)
        prior[u+2] += priorparam->b[n]
               * pow( fabs(x_curr-neighbor[n]), priorparam->p ) ;
    }
    prior[u+2]= prior[u+2]/pow(priorparam->sigma,priorparam->p)/priorparam->p;
    free(neighbor);
  }

  S=physparam->S;
  printf("alpha %e, prior[0] %f,  prior[1] %f, prior[2] %f, prior[3] %f\n",alpha, prior[0], prior[1], prior[2], prior[3]);
  return(S+S*log(alpha)+.5*prior[0]+.5*prior[1]+.5*prior[2]+.5*prior[3]);
}


mtype obj_fn(
  phys_param *physparam,
  prior_param *prior_D,
  prior_param *prior_mua,
  mtype ****x,
  mtype alpha
)
{
  int u,i,j,l,n,S;
  mtype x_curr;
  mtype *neighbor;
  mtype prior[2];
  prior_param *priorparam;
  

  for(u=0; u<2; u++){
    if(u==0){
      priorparam=prior_mua;
    }
    else{
      priorparam=prior_D;
    }

    neighbor = (mtype *) malloc(sizeof(mtype)*priorparam->Nneighbor);

    prior[u]=0.0;
    for(i=1; i<physparam->Ni-1; i++)
    for(j=1; j<physparam->Nj-1; j++)
    for(l=1; l<physparam->Nl-1; l++){
      x_curr=x[u][i][j][l];
      get_neighbor(u,i,j,l, priorparam, x, neighbor);
      for (n=0; n<priorparam->Nneighbor; n++)
        prior[u] += priorparam->b[n]
               * pow( fabs(x_curr-neighbor[n]), priorparam->p ) ;
    }
    prior[u]= prior[u]/pow(priorparam->sigma,priorparam->p)/priorparam->p;
    free(neighbor);
  }

  S=physparam->S;

  return(S+S*log(alpha)+.5*prior[0]+.5*prior[1]);
}



prior_param *alloc_prior_param(
  double (*rho)(double,double,double,double), /* potential fn. rho(xi,xj,sigma,p) */
  dtype sigma, dtype p, 
  int Nneighbor, dtype *b)
{
  int i;
  prior_param *priorparam;
	
  priorparam = (prior_param *) malloc(sizeof(prior_param));
  priorparam->rho   = rho;
  priorparam->sigma = sigma;
  priorparam->p     = p;  
  priorparam->Nneighbor = Nneighbor;

  priorparam->b = (dtype *) malloc(sizeof(dtype)*Nneighbor);  
  for (i=0; i<Nneighbor; i++)
    priorparam->b[i] = b[i];  

  return(priorparam);
}


void free_prior_param(prior_param *priorparam)
{
  free(priorparam->b);
  free(priorparam); 
}


double calc_rmse(
  mtype **yerror,
  phys_param *physparam)
{
  int s;
  double error=0.0;

  for(s=0; s<physparam->S; s++)
    error += AbsSquare(yerror[s][0], yerror[s][1]);

  return(sqrt(error/(physparam->S)));
}


double calc_rmse_with_lambda(
  mtype **yerror,
  mtype *lambda,
  phys_param *physparam)
{
  int s;
  double error=0.0;

  for(s=0; s<physparam->S; s++)
    error += lambda[s] * AbsSquare(yerror[s][0], yerror[s][1]);

  return(sqrt(error/(physparam->S)));
}


void calc_lambda(
  mtype **y,
  mtype **fx, 
  phys_param *physparam,
  src_param *srcparam,
  det_param *detparam,
  mtype *lambda,  /* output */
  mtype *alpha1)    /* output */ 
{
  int s,k,m;
  double alpha;
  double temp1r, temp1i;  
 
  /* compute alpha */
  alpha = 0.0;

  for(s=0; s<physparam->S; s++)
  {
      k = findk(physparam, s);
      m = findm(physparam, s);

      temp1r=calimultr(fx[s][0], fx[s][1], srcparam, detparam, physparam,k,m);
      temp1i=calimulti(fx[s][0], fx[s][1], srcparam, detparam, physparam,k,m);


      alpha += AbsSquare(y[s][0] - temp1r,
                         y[s][1] - temp1i)
              / sqrt(AbsSquare(y[s][0], y[s][1]));
  }
  alpha /= (physparam->S);

  /* compute lambda */
  for(s=0; s<physparam->S; s++) {
    lambda[s] = 1.0 / (alpha * sqrt(AbsSquare(y[s][0], y[s][1]))); 
  }
  *alpha1=alpha;
}      

/* compute theta1 and theta2 for one-pixel ICD update */
/*  1) theta1 = -2 * Re{At^H * lambda * yerror}    */
/*  2) theta2 = 2 * At^H * lambda * At            */	 
void ICD_params_r(
  mtype **At_r,           /* At[s,c]        */
  mtype **yerror_r,      /* yerror[s,c]    */
  mtype *lambda_r,       /* lambda[s]      */
  phys_param *physparam,
  src_param *srcparam,
  det_param *detparam,
  mtype *theta1,       /* output           */
  mtype *theta2,
  int *p_to_s,
  int P,
  int *k_list,
  int *m_list
  )     
{
  int s,p,k,m;
  double temp1r, temp1i, temp2r, temp2i;
  /* Milstein */
  double temp3r, temp3i; 
 
  *theta1 = 0.0;
  *theta2 = 0.0;

  for(p=0; p<P; p++)
  {
      s=p_to_s[p];
      k = k_list[s];
      m = m_list[s];

      temp3r=calimultr(At_r[p][0], At_r[p][1], srcparam, detparam, physparam,k,m);
      temp3i=calimulti(At_r[p][0], At_r[p][1], srcparam, detparam, physparam,k,m);
 
     *theta1 += (temp3r * yerror_r[p][0] + temp3i * yerror_r[p][1]) * lambda_r[p];
     *theta2 += AbsSquare(temp3r, temp3i) * lambda_r[p];
   /* Milstein - End */  
  }

  *theta1 *= (-2.0);
  *theta2 *= (2.0);
}

/* compute theta1 and theta2 for one-pixel ICD update */
/*  1) theta1 = -2 * Re{At^H * lambda * yerror}    */
/*  2) theta2 = 2 * At^H * lambda * At            */	 
void ICD_params(
  mtype **At,           /* At[s,c]        */
  mtype **yerror,      /* yerror[s,c]    */
  mtype *lambda,       /* lambda[s]      */
  phys_param *physparam,
  src_param *srcparam,
  det_param *detparam,
  double *theta1,       /* output           */
  double *theta2)       /* output           */ 
{
  int s,k,m;
  double temp1r, temp1i, temp2r, temp2i;
  /* Milstein */
  double temp3r, temp3i; 
 
  *theta1 = 0.0;
  *theta2 = 0.0;

  for(s=0; s<physparam->S; s++)
  {
      k = findk(physparam, s);
      m = findm(physparam, s);

      temp3r=calimultr(At[s][0], At[s][1], srcparam, detparam, physparam,k,m);
      temp3i=calimulti(At[s][0], At[s][1], srcparam, detparam, physparam,k,m);
 
     *theta1 += (temp3r * yerror[s][0] + temp3i * yerror[s][1]) * lambda[s];
     *theta2 += AbsSquare(temp3r, temp3i) * lambda[s];
   /* Milstein - End */  
  }

  *theta1 *= (-2.0);
  *theta2 *= (2.0);
}
 

int get_neighbor(
  int u, int i, int j, int l,  /* position of pixel being updated */
  prior_param *priorparam,
  mtype ****x,                 /* x[][u,i,j,l]   */
  mtype *neighbor)
{
  int p,q,r,ss=0;

  if (priorparam->Nneighbor==8) {
    for (p=-1; p<2; p++)
    for (q=-1; q<2; q++)
      if (!(p==0 && q==0)) 
        neighbor[ss++] = x[u][i+p][j+q][l]; 
  }
  else if (priorparam->Nneighbor==26) {
    for (p=-1; p<2; p++)
    for (q=-1; q<2; q++)
    for (r=-1; r<2; r++)
      if (!(p==0 && q==0 && r==0)) 
        neighbor[ss++] = x[u][i+p][j+q][l+r]; 
  }
  else {
/*    fprintf(stderr, "Invalid number of neighborhood pixels.\n"); */
    printf("Invalid number of neighborhood pixels.\n");
    exit(1);
  } 

  return 0; 
}


double dev_obj_fn(
  mtype x,                       /* point being evaluated */ 
  mtype x_prev,                  /* value of last iteration */
  prior_param *priorparam,
  double theta1, 
  double theta2, 
  mtype *neighbor)
{
  double prior = 0.0;
  int i; 

  /* calculate prior-model-related part of derivative of object function */ 
  for (i=0; i<priorparam->Nneighbor; i++) 
    prior += priorparam->b[i] 
             * pow( fabs(x-neighbor[i]), priorparam->p-1 ) * SG(x,neighbor[i]); 

  return(theta1+theta2*(x-x_prev)+prior/pow(priorparam->sigma,priorparam->p));
}


/* find root of (derivative of object fn)=0 with a half-interval search */
double  root_find(
  double (*fp)(mtype, mtype, prior_param *, double, double, mtype *), 
  double a,                /* minimum value of solution */
  double b,                /* maximum value of solution */
  double err,              /* accuarcy of solution */
  int *code,               /* error code */
  prior_param *priorparam, /* From this, the parameters are for (*f)() */
  mtype x_prev,           /* value of last iteration */
  double theta1, 
  double theta2, 
  mtype *neighbor)         /* in 3D, neighbor[26] raster-ordered from [-1,-1,-1] to [1,1,1] */
/* Solves equation (*f)(x,constants) = 0 on x in [a,b]. Uses half interval method.*/
/* In ODT, (*f)() will point dev_obj_fn().                                        */        
/* Requires that (*f)(a) and (*f)(b) have opposite signs.                         */
/* Returns code=0  if signs are opposite.                                         */
/* Returns code=1  if signs are both positive.                                    */
/* Returns code=-1 if signs are both positive and code=-1 for both negative.      */
{
  int     signa,signb,signc;
  double  fa,fb,fc,c,  signaling_nan();
  double  dist;

  fa = dev_obj_fn(a, x_prev, priorparam, theta1, theta2, neighbor); 
  signa = fa>0;

  fb = dev_obj_fn(b, x_prev, priorparam, theta1, theta2, neighbor); 
  signb = fb>0;

  /* check starting conditions */
  if( signa==signb ) {
    if(signa==1) *code =  1;
    else         *code = -1;
    return(0.0);
  }
  else *code = 0;


  /* half interval search */
  if( (dist=b-a)<0 )
    dist = -dist;

  while(dist>err) {
    c = (b+a)/2;
    fc = dev_obj_fn(c, x_prev, priorparam, theta1, theta2, neighbor); 
    signc = fc>0;

    if(signa == signc)
    {
      a = c;
      fa = fc;
    }
    else
    {
      b = c;
      fb = fc;
    }

    if( (dist=b-a)<0 )
      dist = -dist;
  }

  /* linear interpolation */
  if( (fb-fa)==0 )
    return(a);
  else {
    c = (a*fb - b*fa)/(fb-fa);
    return(c);
  }
}


double ICD_pixel_update(
  int u, int i, int j, int l,  /* position of pixel being updated */
  mtype **At,                 /* At[u,i,j,l][c,s] for given u,i,j,l  */
  mtype ****x,                 /* x[][u,i,j,l]   */
  mtype **yerror,             /* yerror[s,c]  */
  mtype *lambda,              /* lambda[s]   */
  prior_param *priorparam,
  phys_param *physparam,
  src_param *srcparam,
  det_param *detparam
)
{
  double theta1, theta2;
  double a, b; 
  double xhat;
  mtype  *neighbor; 
  int    n, errcode; 

  neighbor = (mtype *) malloc(sizeof(mtype)*priorparam->Nneighbor);

  /* Compute theta1 and theta2 by calling ICD_params(). */

  ICD_params(At, yerror, lambda, physparam, srcparam, detparam, &theta1, &theta2);

  /* Get neighborhood.  */
  get_neighbor(u,i,j,l, priorparam, x, neighbor);
 
  /* set half interval search range [a,b] for find_root().  */
  if(theta2==0.0) {
    fprintf(stderr, "Divided by theta2=0 in root_df().\n");
    return(x[u][i][j][l]);
  }

  a = x[u][i][j][l] - theta1/theta2;
  b = x[u][i][j][l] - theta1/theta2;

  for(n=0;n<priorparam->Nneighbor;n++)
  {
    if(a > neighbor[n])
      a = neighbor[n];
    else if(b < neighbor[n])
      b = neighbor[n];
  }
  /* root_find().   */
  xhat = root_find(dev_obj_fn, a, b, 0.00001, &errcode,  
                   priorparam, x[u][i][j][l], theta1, theta2, neighbor); 


  /* Some error handling for errcode!=0 or negative xhat */
  if (errcode!=0) {
    xhat = x[u][i][j][l]; 
/*    fprintf(stderr, "error in root_find().\n");   */
    printf("1"); 
  }

  if (xhat <0.0) {
    xhat = 0.0;
/*     printf("DEBUG: xhat less than zero\n");  */

  } 
/*   fprintf(stderr, "mu[%2d,%2d,%2d,%2d] :  %f -> %f,   %e    t1=%e   t2=%e\n",  
          u, i, j, l, x[u][i][j][l], xhat, x[u][i][j][l]-theta1/theta2, 
          theta1, theta2); 
*/   
  free(neighbor);

  return(xhat); 
}

double ICD_pixel_update_fluo_r(
  mtype ****x,                 /* x[][u,i,j,l]   */
  mtype **theta1,              /* theta1[r]   */
  mtype **theta2,              /* theta2[r]   */
  mtype *freqs,
  int u, int i, int j, int l, int R,
  prior_param *priorparam,
  phys_param *physparam,
  src_param *srcparam,
  det_param *detparam,
  mtype *neighbor
)
{

  double a, b; 
  double xhat;
  int    n, errcode; 
  arglist args;

/*  neighbor = (mtype *) malloc(sizeof(mtype)*priorparam->Nneighbor);*/


  /* Get neighborhood.  */
  get_neighbor(u,i,j,l, priorparam, x, neighbor);
 
  /* set half interval search range [a,b] for find_root().  */

    args.priorparam=priorparam;
    args.x=x;
    args.u=u;
    args.i=i;
    args.j=j;
    args.l=l;
    args.R=R;
    args.theta1=theta1;
    args.theta2=theta2;
    args.neighbor=neighbor;
    args.freqs=freqs;


  if(u==0){  
/*      xhat=bracket(0.0, 4e-9, 20, &args, ff)  ;  */
      xhat=bracket(0.0, 4e-7, 30, &args, ff)  ; 
  }
  else{
/*    xhat=bracket(0.0, 0.0010, 40, &args, ff)  ;  */
    xhat=bracket(0.0, 50.000, 55, &args, ff)  ; 
  }


  return(xhat); 
}

double ICD_pixel_update_fluo(
  int u, int i, int j, int l,  /* position of pixel being updated */
  mtype **At,                 /* At[u,i,j,l][c,s] for given u,i,j,l  */
  mtype ****x,                 /* x[][u,i,j,l]   */
  mtype **yerror,             /* yerror[s,c]  */
  mtype *lambda,              /* lambda[s]   */
  prior_param *priorparam,
  phys_param *physparam,
  src_param *srcparam,
  det_param *detparam
)
{
  double theta1, theta2;
  double a, b; 
  double xhat;
  mtype  *neighbor; 
  int    n, errcode; 

  neighbor = (mtype *) malloc(sizeof(mtype)*priorparam->Nneighbor);

  /* Compute theta1 and theta2 by calling ICD_params(). */

  ICD_params(At, yerror, lambda, physparam, srcparam, detparam, &theta1, &theta2);

  /* Get neighborhood.  */
  get_neighbor(u,i,j,l, priorparam, x, neighbor);
 
  /* set half interval search range [a,b] for find_root().  */
  if(fabs(theta2)<1.e-15) {
/*    fprintf(stderr, "Divided by theta2=0 in root_df().\n");
    printf("u=%d\n",u); */
    return(x[u][i][j][l]);
  }

  a = x[u][i][j][l] - theta1/theta2;
  b = x[u][i][j][l] - theta1/theta2;


  for(n=0;n<priorparam->Nneighbor;n++)
  {
    if(a > neighbor[n])
      a = neighbor[n];
    else if(b < neighbor[n])
      b = neighbor[n];
  }
  /* root_find().   */

  if(u==0){  
    xhat = root_find(dev_obj_fn, a, b, 1.0e-9, &errcode,  
                   priorparam, x[u][i][j][l], theta1, theta2, neighbor); 
  }
  else{
    xhat = root_find(dev_obj_fn, a, b, 1.0e-9, &errcode,  
                   priorparam, x[u][i][j][l], theta1, theta2, neighbor); 
  }


  /* Some error handling for errcode!=0 or negative xhat */
  if (errcode!=0) {
    xhat = x[u][i][j][l]; 
/*    fprintf(stderr, "error in root_find().\n");   */
    printf("1"); 
  }

  if (xhat <0.0) {
    xhat = 0.0;
  } 

  free(neighbor);

  return(xhat); 
}


void ICD_update(
  mtype ****x,          /* initial value x[u][i][j][l] and output */
  mtype **yerror,      /* initial error and output               */   
  mtype *****phi,         /* phi[k][i][j][l][c]                     */ 
  mtype *****green,       /* green[m][i][j][l][c]                   */
  /* BEGIN-NEWCALI */
  mtype **y,
  mtype **fx,
  /* END-NEWCALI */
  mtype *lambda,        
  prior_param *priorparam_D, 
  prior_param *priorparam_mua, 
  src_param *srcparam,
  det_param *detparam,
  phys_param *physparam,
  config_param *config)
{
  int u,i,j,l,ubegin,uend;
  int q,t,rand_index;
  int Ni,Nj,Nl;
  int Nq;
  mtype **At;
  double xhat;
  int s;
  int *rand_update_mask=NULL;
  double temp1r, temp1i, temp2r, temp2i;
  int k,m;
  double temp3r, temp3i;    

  /* BEGIN-NEWCALI */
  double mean[2];
  /* END-NEWCALI */  

  Ni=physparam->Ni;
  Nj=physparam->Nj;
  Nl=physparam->Nl;
  Nq=(Ni-2*config->borderi)*(Nj-2*config->borderj)*(Nl-2*config->borderl);
  
  At = multialloc(sizeof(mtype), 2, physparam->S, 2);

  /* BEGIN-NEWCALI *
  fx = multialloc(sizeof(double), 2, physparam->S, 2);
 
  for (s=0; s<physparam->S; s++)
  for (c=0; c<2; c++) 
    fx[s][c] = y[s][c] - yerror[s][c];
  
  * END-NEWCALI */
                
  rand_update_mask = (int *)malloc(Nq*sizeof(int));

  if (config->mua_flag && config->D_flag) {
    ubegin = 0; 
    uend = 1;
  } 
  else if (config->mua_flag && !config->D_flag) {
    ubegin = 0; 
    uend = 0;
  } 
  else if (!config->mua_flag && config->D_flag) {
    ubegin = 1; 
    uend = 1;
  }
  else {
/*    fprintf(stderr,"Error! Both mua_flag and D_flag are 0!!!\n");
    exit(1); */
    free(rand_update_mask);
    multifree(At, 2);
/*    return 0; */
  }


  for (u=ubegin; u<uend+1; u++){     
    for (q=0; q<Nq; q++){
      rand_update_mask[q]=0;
    }
    for (t=0; t<Nq; t++){
      srand(time(NULL));
      rand_index=rand()%(Nq-t);
      q=-1;
      while(rand_index>=0){
        q++;
        while(rand_update_mask[q]) q++;
        rand_index--;
      }
      rand_update_mask[q]=1;

      l = config->borderl + q%(Nl-2*config->borderl);
      j = config->borderj + (q/(Nl-2*config->borderl)) % (Nj-2*config->borderj);
      i = config->borderi + (q/((Nl-2*config->borderl) * (Nj-2*config->borderj)))
          % (Ni-2*config->borderi);

      calc_frechet_col(u, i, j, l, phi, green, x,
                     srcparam, physparam, At);
      if(u==MUA)
        xhat = ICD_pixel_update(u, i, j, l, At, x, yerror, lambda,
                                priorparam_mua, physparam, srcparam, detparam);
      else
        xhat = ICD_pixel_update(u, i, j, l, At, x, yerror, lambda,
                                priorparam_D, physparam, srcparam, detparam);

      /* update yerror: yerror -= At[u,i,j,l] * delta_x.   */

      if (config->homogeneous_flag==0) {
        for (s=0; s<physparam->S; s++) {
          k = findk(physparam, s);
          m = findm(physparam, s);

          temp3r=calimultr(At[s][0], At[s][1], srcparam, detparam, physparam,k,m);
          temp3i=calimulti(At[s][0], At[s][1], srcparam, detparam, physparam,k,m);
   
          yerror[s][0] -= temp3r * (xhat-x[u][i][j][l]);
          yerror[s][1] -= temp3i * (xhat-x[u][i][j][l]);
        }
      }

      /* update pixel */
      x[u][i][j][l] = xhat;
    }
    printf("yerror=(%e %e)    At=%e\n", yerror[0][0],yerror[0][1], At[0][0]); 
  /* Crude debugging */
    for (q=0; q<Nq; q++){
      if(rand_update_mask[q]!=1){
         printf("Error!!!!!!!!!!!!!!!!!!\n");
         exit(1);
      }
    }
  /*End Crude debugging*/
  }

  /* BEGIN-TEST : Fill all voxels with mean value */
  if (config->homogeneous_flag==1) {
    mean[0] = 0.0;
    mean[1] = 0.0;

    for (u=0; u<2; u++)
    for (i= config->borderi; i<Ni-config->borderi; i++)
    for (j= config->borderj; j<Nj-config->borderj; j++)
    for (l= config->borderl; l<Nl-config->borderl; l++) {
      mean[u] += x[u][i][j][l];
    }

    for (u=0; u<2; u++) {
      mean[u] /= (Nl-2*config->borderl) * (Nj-2*config->borderj) * (Ni-2*config->borderi);

      for (i= 0; i<Ni; i++)
      for (j= 0; j<Nj; j++)
      for (l= 0; l<Nl; l++) {
        for (s=0; s<physparam->S; s++) {
          k = findk(physparam, s);
          m = findm(physparam, s);
          
	  temp3r=calimultr(At[s][0], At[s][1], srcparam, detparam, physparam,k,m);
          temp3i=calimulti(At[s][0], At[s][1], srcparam, detparam, physparam,k,m);
          

          yerror[s][0] -= temp3r * (mean[u]-x[u][i][j][l]);
          yerror[s][1] -= temp3i * (mean[u]-x[u][i][j][l]);
        }
      }
    }

    fprintf(stderr, "%f    %f \n", mean[0], mean[1]); 
  }
  /* END-TEST */

  free(rand_update_mask);
  multifree(At, 2);
}                              

/* Reverse solver using single-grid ICD algorithm.
   This will be used in multigrid reverse solver subroutine rvssolver_mg().
   Output : x (that is, mu_a & mu_s) and yerror  */
void rvssolver_sg(
  mtype ****x,                 /* initial value of mu_a & mu_s and output */
  mtype **y,                  /* observed phi, that is, desired value of f(x^) */
  phys_param *physparam,       /* physics parameters */
  src_param  *srcparam,        /* src_param array */
  det_param  *detparam,        /* detpos array    */
  prior_param *priorparam_D,
  prior_param *priorparam_mua,
  mtype **yerror,             /* output : observ - f(x^)  */
  config_param *config)
{
  int iter, s, c;
  mtype **fx, **fx2;        /* Calculated phi at measurement positions: fx[c,s]*/
  mtype *lambda;
  mtype alpha,cost;
  mtype *****phi;
  mtype *****green;
  double yrmse=9999999.0;
  char filename[255],arrayname[100];
  FILE *fp;
  FILE *fpresult;
  double temp1r, temp1i;
  int k,m;
  /* Milstein */
  double temp2r, temp2i;
  double oldwgtr, oldwgti;
 
  fx     = multialloc(sizeof(mtype), 2, physparam->S, 2);
  lambda = (mtype *)malloc((physparam->S)*sizeof(mtype));
 
  phi    = multialloc(sizeof(mtype), 5, physparam->K,
                      physparam->Ni, physparam->Nj, physparam->Nl, 2);
  green  = multialloc(sizeof(mtype), 5, physparam->M,
                      physparam->Ni, physparam->Nj, physparam->Nl, 2);
 
  for (iter=0; iter<config->niterations && yrmse>config->rmse_tol; iter++)
  {
    printf("CALC_PHI...\n");
    calc_phi(srcparam, detparam, x, physparam, phi, fx);
    printf("calc_green...\n");
    calc_green(srcparam, detparam, x, physparam, green);
 
    sprintf(filename, "%s/RESULT",config->resultpath);
    fpresult = fopen(filename, "a");
    fprintf(fpresult, "%2d \t", iter+1);
 
    calc_lambda(y, fx, physparam, srcparam, detparam, lambda, &alpha);
 
    cost=obj_fn( physparam, priorparam_D, priorparam_mua, x, alpha);
 
    fprintf(fpresult, "%e  \t%e  \t", cost, alpha);
    fprintf(fpresult, "%e  \t", calc_rmse(yerror, physparam));
 
    for (s=0; s<physparam->S; s++) {
      k = findk(physparam, s);
      m = findm(physparam, s);

      temp1r=calimultr(fx[s][0], fx[s][1], srcparam, detparam, physparam,k,m);
      temp1i=calimulti(fx[s][0], fx[s][1], srcparam, detparam, physparam,k,m);
      
      yerror[s][0] = y[s][0] - temp1r;
      yerror[s][1] = y[s][1] - temp1i;        
    }
 
    if (config->calibration_flag == 1)
      src_det_calibrate(y, fx, lambda, srcparam, detparam, physparam);
 
    if (config->global_weight_flag == 1)
      global_weight_calibrate(x, y, fx, lambda, srcparam, detparam, physparam);
   
    if (config->calibration_flag == 1 || config->global_weight_flag == 1) {
      for (s=0; s<physparam->S; s++) {
        k = findk(physparam, s);
        m = findm(physparam, s);

        temp1r=calimultr(fx[s][0], fx[s][1], srcparam, detparam, physparam,k,m);
        temp1i=calimulti(fx[s][0], fx[s][1], srcparam, detparam, physparam,k,m);
      
        yerror[s][0] = y[s][0] - temp1r;
        yerror[s][1] = y[s][1] - temp1i;        

      }
    }

    ICD_update(x, yerror, phi, green, y, fx,
               lambda, priorparam_D, priorparam_mua,
               srcparam, detparam, physparam, config);

/*                     
    Homo_update_Newton(x, yerror, phi, green, lambda, srcparam, detparam, physparam, config);
    Homo_update_Search2(x, y, fx, lambda, phi, green, yerror, srcparam, detparam, physparam, config);
*/
 
    yrmse = calc_rmse(yerror, physparam);
 
    fprintf(fpresult, "%e  \n", yrmse);
    fclose(fpresult);
 
    if (config->mu_store_flag !=0) {
      sprintf(filename, "%s/muhat%d.dat", config->muhatpath, iter);
      sprintf(arrayname, "muhat%d", iter);
      fp = datOpen(filename, "w+b");
      write_float_array(fp, arrayname, &x[0][0][0][0], 4, 2,
                        physparam->Ni, physparam->Nj, physparam->Nl);
      datClose(fp);
    }
 
    sprintf(filename, "%s/WEIGHTS",config->resultpath);
    fpresult = fopen(filename, "a");
    fprintf(fpresult, "iter : %2d \n", iter+1);

    if (config->global_weight_flag == 1)
      fprintf(fpresult, "global weight = %f + %f j        D = %f \n",  
              physparam->wgtr, physparam->wgti, x[1][8][8][8]);

    if (config->calibration_flag == 1) {
      for (k=0; k<physparam->K;k++)
        fprintf(fpresult, "source %2d :   %f  %f \n",  k, srcparam[k].calir, srcparam[k].calii);
 
      for (m=0; m<physparam->M;m++)
        fprintf(fpresult, "detector %2d :  %f  %f \n", m, detparam[m].calir, detparam[m].calii);
    }
    fclose(fpresult);

  }
 
  multifree(fx, 2);
  free(lambda);
  multifree(phi, 5);
  multifree(green, 5);
}                        


/* Reverse solver using multigrid will be considered in the future.
   This calls rvssolver_sg() */
void rvssolver_mg();

int findr(src_param *src, int karg){

   double *freqs;
   int numfreqs=0, k, r, found_flag=0;

   freqs=(double *)malloc(karg*sizeof(double));

   for(k=0; k<=karg; k++){
     found_flag=0;
     for(r=0; r<numfreqs; r++){
       if(fabs(freqs[r]-src[k].omega)<(1.e-4*freqs[r])){
         found_flag=1;
	 break;
       }

     }
     if(!found_flag){
	freqs[numfreqs]=src[k].omega;
        numfreqs++;
     }
   }
   free(freqs);
   return (numfreqs-1);
}

void ICD_update_fluo(
  mtype ****fluo,          /* initial value fluo[u][i][j][l] and output */
  mtype **yerror,      /* initial error and output               */   
  mtype *****phi_x,       /* phi_x[k][i][j][l][c]                     */ 
  mtype *****green,       /* green[m][i][j][l][c]                   */
  /* BEGIN-NEWCALI */
  mtype **y,
  mtype **meas_m,
  /* END-NEWCALI */
  mtype *lambda,        
  prior_param *priorparam_tau, 
  prior_param *priorparam_gamma, 
  src_param *srcparam,
  det_param *detparam,
  phys_param *physparam,
  fluo_config_param *config)
{
  int u,i,j,l,ubegin,uend;
  int q,t,rand_index;
  int Ni,Nj,Nl;
  int Nq;
  mtype ***At, **At_r, **yerror_r, *lambda_r, *freqs;
  mtype **theta1, **theta2;
  int *p_to_s, *P, R, u1, r, p;
  int *r_list, *k_list, *m_list;
  double fluohat, fluo_lin_r, fluo_lin_i, tau0, tau1, eta0, eta1;
  double omega; 
  int s;
  int *rand_update_mask=NULL;
  double temp1r, temp1i, temp2r, temp2i;
  int k,m;
  double temp3r, temp3i;    
  mtype *neighbor;

  Ni=physparam->Ni;
  Nj=physparam->Nj;
  Nl=physparam->Nl;
  Nq=(Ni-2*config->borderi)*(Nj-2*config->borderj)*(Nl-2*config->borderl);
 
  At = multialloc(sizeof(mtype), 3, 2, physparam->S, 2);
  At_r = multialloc(sizeof(mtype), 2, physparam->S, 2);
  yerror_r = multialloc(sizeof(mtype), 2, physparam->S, 2);
  lambda_r = (mtype *)malloc((physparam->S)*sizeof(mtype));
  p_to_s=(int *)malloc((physparam->S)*sizeof(int));
  r_list=(int *)malloc((physparam->S)*sizeof(int));
  k_list=(int *)malloc((physparam->S)*sizeof(int));
  m_list=(int *)malloc((physparam->S)*sizeof(int));
  P=(int *)malloc((physparam->S)*sizeof(int));

  R=findr(srcparam,(physparam->K)-1)+1;
  printf("R=%d\n",R);
  freqs=(mtype *)malloc(R*sizeof(mtype));
  theta1=multialloc(sizeof(mtype), 2, 2, R);
  theta2=multialloc(sizeof(mtype), 2, 2, R);

  rand_update_mask = (int *)malloc(Nq*sizeof(int));
  
  neighbor = (mtype *) malloc(sizeof(mtype)*priorparam_gamma->Nneighbor);

  if (config->tau_flag && config->etamu_flag) {
    ubegin = 0; 
    uend = 1;
  } 
  else if (config->tau_flag && !config->etamu_flag) {
    ubegin = 0; 
    uend = 0;
  } 
  else if (!config->tau_flag && config->etamu_flag) {
    ubegin = 1; 
    uend = 1;
  }
  else {
/*    fprintf(stderr,"Error! Both tau_flag and etamu_flag are 0!!!\n"); */
     printf("Error! Both tau_flag and etamu_flag are 0!!!\n"); 
    exit(1);
  }

  for(s=0; s<physparam->S; s++){
     k_list[s]=findk(physparam,s); 
     m_list[s]=findm(physparam,s); 
     r_list[s]=findr(srcparam,k_list[s]); 
  }

  for (u=ubegin; u<uend+1; u++){     
    for (q=0; q<Nq; q++){
      rand_update_mask[q]=0;
    }
    for (t=0; t<Nq; t++){
      srand(time(NULL));
      rand_index=rand()%(Nq-t);
      q=-1;
      while(rand_index>=0){
        q++;
        while(rand_update_mask[q]) q++;
        rand_index--;
      }
      rand_update_mask[q]=1;

      l = config->borderl + q%(Nl-2*config->borderl);
      j = config->borderj + (q/(Nl-2*config->borderl)) % (Nj-2*config->borderj);
      i = config->borderi + (q/((Nl-2*config->borderl) * (Nj-2*config->borderj)))
        % (Ni-2*config->borderi);

      for(u1=0; u1<=1; u1++){

        calc_frechet_col_fluo(u1, i, j, l, phi_x, green, fluo,
                   srcparam, physparam, At[u1]);

        for(r=0; r<R; r++){
          p=0; 
          for(s=0; s<physparam->S; s++){
            if(r_list[s]==r){
              freqs[r]=srcparam[k_list[s]].omega;
              At_r[p][0]=At[u1][s][0];
              At_r[p][1]=At[u1][s][1];
              yerror_r[p][0]=yerror[s][0];
              yerror_r[p][1]=yerror[s][1];
              lambda_r[p]=lambda[s];
	      p_to_s[p]=s;
	      p++;
            }
          }
          P[r]=p; 
	  ICD_params_r(At_r, yerror_r, lambda_r, physparam, srcparam, detparam, 
			   &theta1[u1][r], &theta2[u1][r],p_to_s, p, k_list, m_list); 

        }
      }

      if(u==0){
         fluohat=ICD_pixel_update_fluo_r(
           fluo, theta1, theta2, freqs,   u, i, j, l,  R,
           priorparam_tau, physparam, srcparam, detparam,neighbor);

      }
       else{
         fluohat=ICD_pixel_update_fluo_r(
           fluo, theta1, theta2, freqs,   u, i, j, l,  R,
           priorparam_gamma, physparam, srcparam, detparam,neighbor);
       }

      if (config->homogeneous_flag==0) {
        for (s=0; s<physparam->S; s++) {
          k = findk(physparam, s);
          m = findm(physparam, s);
	  omega=srcparam[k].omega;
	  tau0=fluo[0][i][j][l];
	  eta0=fluo[1][i][j][l];
	  
	  if(u==0){
	    tau1=fluohat;
	    eta1=eta0;
          }
	  else{
	    tau1=tau0;
	    eta1=fluohat;
          }
          
       	  fluo_lin_r= ( eta1/(1.+omega*omega*tau1*tau1) ) - ( eta0/(1.+omega*omega*tau0*tau0) );
          fluo_lin_i= ( -omega*tau1*eta1/(1.+omega*omega*tau1*tau1) ) - ( -omega*tau0*eta0/(1.+omega*omega*tau0*tau0) );

          temp3r=calimultr(At[0][s][0], At[0][s][1], srcparam, detparam, physparam,k,m);
          temp3i=calimulti(At[0][s][0], At[0][s][1], srcparam, detparam, physparam,k,m);

          yerror[s][0] -= cplxmultr(temp3r, temp3i, fluo_lin_r, fluo_lin_i);
          yerror[s][1] -= cplxmulti(temp3r, temp3i, fluo_lin_r, fluo_lin_i);
          meas_m[s][0] += cplxmultr(temp3r, temp3i, fluo_lin_r, fluo_lin_i);
          meas_m[s][1] += cplxmulti(temp3r, temp3i, fluo_lin_r, fluo_lin_i);
 
        }
      }
      /* update pixel */
      fluo[u][i][j][l] = fluohat;
    }
    printf("yerror=(%e %e)    At=%e\n", yerror[0][0],yerror[0][1], At[0][0][0]); 
  /* Crude debugging */
    for (q=0; q<Nq; q++){
      if(rand_update_mask[q]!=1){
         printf("Error!!!!!!!!!!!!!!!!!!\n");
         exit(1);
      }
    }
  /*End Crude debugging*/
  }

  free(rand_update_mask);
  free(freqs);
  free(p_to_s);
  free(k_list);
  free(m_list);
  free(r_list);
  free(neighbor);
  free(P);
  multifree(theta1,2);
  multifree(theta2,2);
  free(lambda_r);
  multifree(yerror_r,2);
  multifree(At, 3);
  multifree(At_r, 2);
}                              


void rvssolver_sg_fluo(
  mtype ****fluo,                 /* initial value of tau and etamu and output */\
   mtype ****mu_m,
   mtype *****phi_x,
   mtype *****phi_m,
   mtype **meas_m,
  mtype **y,                  /* observed phi, that is, desired value of f(x^) */
  phys_param *physparam,       /* physics parameters */
  src_param  *srcparam,        /* src_param array */
  det_param  *detparam,        /* detpos array    */
  prior_param *priorparam_tau,
  prior_param *priorparam_gamma,
  prior_param *priorparam_D,
  prior_param *priorparam_mua,
  mtype **yerror,             /* output : observ - f(x^)  */
  fluo_config_param *config)
{
  int iter, s;
  mtype *lambda;
  mtype alpha,cost;
  mtype *****green;
  double yrmse=9999999.0;
  char filename[255],arrayname[100];
  FILE *fp;
  FILE *fpresult;
  double temp1r, temp1i;
  int k,m;
  /* Milstein */
  double temp2r, temp2i;
  config_param config_x;
  int fluo_updates=0;

  convert_config(config, &config_x);
 
  lambda = (mtype *)malloc((physparam->S)*sizeof(mtype));
 
  green  = multialloc(sizeof(mtype), 5, physparam->M,
                      physparam->Ni, physparam->Nj, physparam->Nl, 2);

  printf("calc_green...\n");
  calc_green(srcparam, detparam, mu_m, physparam, green);

  for (iter=0; iter<config->niterations && yrmse>config->rmse_tol; iter++)
/*  for (iter=0; iter<1; iter++) */
  {
    if(iter!=0 && (config->mua_flag || config->D_flag) ){
      printf("calc_phi...\n");
      calc_phi_fluo(fluo, phi_x, srcparam, detparam, mu_m, physparam, phi_m, meas_m); 
      printf("calc_green...\n");
      calc_green(srcparam, detparam, mu_m, physparam, green); 
    }
/*
    mult_frechet_row(phi_x,green,fluo, srcparam,detparam, physparam);
 
    fp = datOpen("phi_x.dat", "w+b");
    write_float_array(fp, "phi_x", &phi_x[0][0][0][0][0], 5, physparam->K, physparam->Ni,
                           physparam->Nj, physparam->Nl, 2);
    datClose(fp);
    fp = datOpen("green.dat", "w+b");
    write_float_array(fp, "green", &green[0][0][0][0][0], 5, physparam->M, physparam->Ni,
                           physparam->Nj, physparam->Nl, 2);
    datClose(fp);

*/
    sprintf(filename, "%s/RESULT",config->resultpath);
    fpresult = fopen(filename, "a");
    fprintf(fpresult, "%2d \t", iter+1);
  
    calc_lambda(y, meas_m, physparam, srcparam, detparam, lambda, &alpha);
 
    cost=obj_fn_fluo( physparam, 
                       priorparam_tau, priorparam_gamma, 
                       priorparam_D , priorparam_mua, 
                       fluo , mu_m,    alpha);
 
    fprintf(fpresult, "%e  \t%e  \n", cost, alpha);
/*    fprintf(fpresult, "%e  \t|", calc_rmse(yerror, physparam)); */

 
    for (s=0; s<physparam->S; s++) {
      k = findk(physparam, s);
      m = findm(physparam, s);
     
      temp1r=calimultr(meas_m[s][0], meas_m[s][1], srcparam, detparam, physparam,k,m);
      temp1i=calimulti(meas_m[s][0], meas_m[s][1], srcparam, detparam, physparam,k,m);
      yerror[s][0] = y[s][0] - temp1r;
      yerror[s][1] = y[s][1] - temp1i;         
    }
 
    if (config->calibration_flag == 1)
      src_det_calibrate(y, meas_m, lambda, srcparam, detparam, physparam);
 
    if (config->global_weight_flag == 1)
      global_weight_calibrate(fluo, y, meas_m, lambda, srcparam, detparam, physparam);
    
    if (config->calibration_flag == 1 || config->global_weight_flag == 1) {
      for (s=0; s<physparam->S; s++) {
        k = findk(physparam, s);
        m = findm(physparam, s);

        temp1r=calimultr(meas_m[s][0], meas_m[s][1], srcparam, detparam, physparam,k,m);
        temp1i=calimulti(meas_m[s][0], meas_m[s][1], srcparam, detparam, physparam,k,m);
        yerror[s][0] = y[s][0] - temp1r;
        yerror[s][1] = y[s][1] - temp1i;
      }
    }

      ICD_update_fluo(fluo, yerror, phi_x, 
                    green, y, meas_m, lambda, priorparam_tau, priorparam_gamma,
                    srcparam, detparam, physparam, config);


    calc_lambda(y, meas_m, physparam, srcparam, detparam, lambda, &alpha);
    cost=obj_fn_fluo( physparam, 
                       priorparam_tau, priorparam_gamma, 
                       priorparam_D , priorparam_mua, 
                       fluo, mu_m,    alpha);
/*    fprintf(fpresult, "AF: %e  \t%e  \t", cost, alpha);
    fprintf(fpresult, "%e  \t|", calc_rmse(yerror, physparam)); */

    yrmse = calc_rmse(yerror, physparam);
 
/*    fprintf(fpresult, "%e  \n", yrmse); */
    fclose(fpresult);
  
    if (config->fluo_store_flag !=0) {
      sprintf(filename, "%s/fluohat.dat", config->fluohatpath);
      sprintf(arrayname, "fluohat");
      fp = datOpen(filename, "w+b");
      write_float_array(fp, arrayname, &fluo[0][0][0][0], 4, 2,
                        physparam->Ni, physparam->Nj, physparam->Nl);
      datClose(fp);
    }
 
    sprintf(filename, "%s/WEIGHTS",config->resultpath);
    fpresult = fopen(filename, "a");
    fprintf(fpresult, "iter : %2d \n", iter+1);


    if (config->calibration_flag == 1) {
      for (k=0; k<physparam->K;k++)
        fprintf(fpresult, "source %2d :   %f  %f \n",  k, srcparam[k].calir, srcparam[k].calii);
 
      for (m=0; m<physparam->M;m++)
        fprintf(fpresult, "detector %2d :  %f  %f \n", m, detparam[m].calir, detparam[m].calii);
    }
    fclose(fpresult);
  }
   
  free(lambda);
  multifree(green, 5);
}                        

int convert_fluo(
   int dir,
   phys_param *physparam,  /* physics parameters */
   src_param  *srcparam,   /* src_param array */
   mtype ****fluo_in,
   mtype ****fluo_out
) 
{
  double gamma, tau, etamu, omega;
  int i, j, l, Ni, Nj, Nl;
  char filename[255],arrayname[100];
  FILE *fp;

  omega=srcparam[0].omega;
  Ni=physparam->Ni;
  Nj=physparam->Nj;
  Nl=physparam->Nl;

  if(dir==1){  /* Convert from tau, etamu to alpha, beta */
    for (i= 0; i<Ni; i++)
    for (j= 0; j<Nj; j++)
    for (l= 0; l<Nl; l++) {
      tau=fluo_in[0][i][j][l];
      etamu=fluo_in[1][i][j][l];
      gamma=etamu/(1.0+omega*omega*tau*tau);
      fluo_out[0][i][j][l]=tau;
      fluo_out[1][i][j][l]=gamma;
    }
/*    sprintf(filename, "%s/alpha.dat", "DUMPFLUO1");
    sprintf(arrayname, "alpha");
    fp = datOpen(filename, "w+b");
    write_float_array(fp, arrayname, &fluo_out[0][0][0][0], 4, 2,
                        physparam->Ni, physparam->Nj, physparam->Nl);
    datClose(fp); */
  }
    
  else{
    for (i= 0; i<Ni; i++)
    for (j= 0; j<Nj; j++)
    for (l= 0; l<Nl; l++) {
      tau=fluo_in[0][i][j][l];
      gamma=fluo_in[1][i][j][l];
   
      etamu=gamma*(1.0+omega*omega*tau*tau);
      fluo_out[0][i][j][l]=tau;
      fluo_out[1][i][j][l]=etamu;
    }
  }
  return 0;
}


int convert_config(fluo_config_param *config, config_param *config_x){

  config_x->borderi= config->borderi;
  config_x->borderj= config->borderj;
  config_x->borderl= config->borderl;

  config_x->niterations = config->niterations;
  config_x->rmse_tol = config-> rmse_tol ;
  config_x-> alpha_bound = config-> alpha_bound ;

  config_x-> mua_backg = config-> mua_backg ;
  config_x-> D_backg   = config-> D_backg ;

  config_x-> init_wgtr = config-> init_wgtr ;
  config_x-> init_wgti = config-> init_wgti ;


  strcpy(config_x-> muhatpath,  config-> fluohatpath);
  config_x-> mu_store_flag = config-> fluo_store_flag ;
  strcpy(config_x-> resultpath,  config-> resultpath);

  strcpy(config_x-> init_guess_path,  config-> init_guess_opt_path);
  strcpy(config_x-> init_guess_varname,  config-> init_guess_opt_varname);

  strcpy(config_x-> phys_file,  config-> phys_file);
  strcpy(config_x-> prior_D_file,  config-> prior_D_file);
  strcpy(config_x-> prior_mua_file,  config-> prior_mua_file);
  strcpy(config_x-> meas_file,  config-> meas_file);

  config_x->calibration_flag   =  config-> calibration_flag  ;
  config_x->homogeneous_flag   =  config-> homogeneous_flag ;
  config_x->mua_flag  =  config->mua_flag  ;
  config_x->D_flag  =  config->D_flag  ;

  config_x-> global_weight_flag =  config-> global_weight_flag ;
  config_x-> init_guess_flag =  config-> init_guess_opt_flag ;
  
  return 0;
}
