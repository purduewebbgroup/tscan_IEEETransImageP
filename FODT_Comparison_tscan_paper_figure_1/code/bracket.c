#include "defs.h"
#include <math.h>
#define min(X, Y)  ((X) < (Y) ? (X) : (Y))
#define max(X, Y)  ((X) > (Y) ? (X) : (Y))

mtype  bracket(
   mtype as,
   mtype bs,
   int N, 
   arglist *args,
   mtype (*f)(arglist *, mtype)
)	/* ANSI: float (*f)(float); */
{
  int k;
  mtype del,fac, feval, fmin=1e12, dexmin, dex , ax, bx, cx, tol;
  mtype xmin;
  tol=1.0e-5;

  /* del=pow(10,log10(bs/as)/N); */
  del=(bs-as)/((double)N); 

/*  fac=1.0; */
  fac=0.0;
  for(k=0; k<=N; k++){
    dex=as+fac;
    feval=(*f)(args,dex);
/*    printf("%f\n",feval); */
    if(feval<fmin){
       dexmin=dex;
       fmin=feval; 
    }
    fac=fac+del;
  }
  ax=max( dexmin-del, as);
  bx=dexmin;
  cx=min( dexmin+del, bs);


  if(ax==bx ){
    golden(ax,ax+del,ax+2.0*del,f,args,tol,36, &xmin);
  }
  else if(cx==bx){
    golden(cx-2.0*del,cx-del,cx,f,args,tol,36, &xmin);
  }
  else{
    golden(ax,bx,cx,f,args,tol,36, &xmin);
    return xmin;
  }

  return xmin;
}

mtype ff( arglist *args, mtype fluohat)
{

  prior_param *priorparam=args->priorparam;
  mtype ****fluo=args->x;
  int u=args->u;
  int i=args->i;
  int j=args->j;
  int l=args->l;
  int R=args->R;
  mtype **theta1=args->theta1;
  mtype **theta2=args->theta2;
  mtype *neighbor=args->neighbor;
  mtype *freqs=args->freqs;

  mtype data=0.0, prior = 0.0, omega;
  mtype tau1, eta1, tau0, eta0, del_fluo_r, del_fluo_i;
  int n,r;

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

  for (r=0; r<R; r++) {
     omega=freqs[r];
     del_fluo_r= ( eta1/(1.+omega*omega*tau1*tau1) ) 
               - ( eta0/(1.+omega*omega*tau0*tau0) );
     del_fluo_i= ( -omega*tau1*eta1/(1.+omega*omega*tau1*tau1) ) 
               - ( -omega*tau0*eta0/(1.+omega*omega*tau0*tau0) );
   
     if(theta1[1][r]*theta2[1][r]<0.0){
/*       printf("theta1[0][r] %f ; theta2[0][r] %f\n", theta1[0][r], theta2[0][r]);
       printf("theta1[1][r] %f ; theta2[1][r] %f\n", theta1[1][r], theta2[1][r]);
       printf("i=%d, j=%d, l=%d\n",i,j,l);
       exit(1);  */
     } 

     data+= theta1[0][r]*del_fluo_r+0.5*theta2[0][r]*del_fluo_r*del_fluo_r
           +theta1[1][r]*del_fluo_i+0.5*theta2[1][r]*del_fluo_i*del_fluo_i;
  }

  for (n=0; n<priorparam->Nneighbor; n++)
       prior += priorparam->b[n] 
	            * pow( fabs(fluohat-neighbor[n]), priorparam->p );
  return (prior/(priorparam->p*pow(priorparam->sigma,priorparam->p))  +data);
  
}
