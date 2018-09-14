/***********<< defs.h >>*********************************************************/

#include "structs.h"
#include "structs2.h"
#include "listfun.h"
#include <stdio.h>
#include <math.h>
#include <stdarg.h>
#include "allocate.h"
#include "randlib.h"
#include "fileop.h"
#include "textfile.h"
#include <time.h>
#include <string.h>

#define  DEFAULT_CYCLES  ( 16 )


/* Constants for index u */
#define MUA (0)
#define MUS (1)

#define  round(x)       ((int)(x+0.5))
#define  SG(A,B)        ((A) > (B) ? 1.0 : -1.0)
#define  AbsSquare(R,I) ((R)*(R)+(I)*(I))
#define  cplxmultr(ar,ai,br,bi)  ((ar)*(br)-(ai)*(bi))
#define  cplxmulti(ar,ai,br,bi)  ((ar)*(bi)+(ai)*(br))



/******************************************************************/
/***                     fwdsolver3d.f                          ***/
/******************************************************************/

/* FORTRAN function declarations  designed for SUN compilation*/
/* Solves elliptic PDE in form of diffusion equation */
extern int fwd3df_(
     /* Input */
  int *nnx,   /* 3-D array dimensions */
  int *nny, 
  int *nnz,
  mtype *xmin, /* Ranges of physical coords */
  mtype *xmax, 
  mtype *ymin, 
  mtype *ymax, 
  mtype *zmin, 
  mtype *zmax, 
  mtype *alpha, /* Array containing 1st PDE coeff: In C, alpha[i][j][l][c] 
                 * In FORTRAN, alpha(c+1,l+1,j+1,i+1) */
  mtype *beta,  /* Array containing 2nd PDE coeff: beta[i][j][l][c] 
                 * In FORTRAN, beta(c+1,l+1,j+1,i+1) */
  mtype *rhs,   /* Array containing RHS of PDE : rhs[i][j][l][c] 
                 * In FORTRAN, rhs(c+1,l+1,j+1,i+1) */
  int *iguessflag, /* 1 if initial guess should be used, 0 otherwise  */
  mtype *mudpkworkspace,
  int *mudpkworksize,
  mtype *rhsworkspace,
  double *phiworkspace,
  int *ncycles,    /* Number of multigrid cycles to be used by MUDPACK
                    * See cud3.d for more details */
  
    /* Input/Output */
  mtype *phi       /* Input: initial guess, if iguessflag=1.  Should be 
                    *   initialized to zeros otherwise.
                    * Output: solution of PDE.  In C, phi[i][j][l][c]
                    * In FORTRAN, phi(c+1,l+1,j+1,i+1) */
);

/* End FORTRAN declarations */

/******************************************************************/
/***                     fwdsolver.c                            ***/
/******************************************************************/
/*
int add_detector_noise(
  phys_param *physparam, 
  dtype alpha_fixed, 
  mtype **meas, 
  mtype *snr);
*/
int add_detector_noise(phys_param *physparam, dtype alpha_fixed,
                                              dtype max_snr, /* Milstein */
                                              mtype **meas, mtype *snr);
int mult_frechet_row(
  /* Input */
  mtype *****phi,
  mtype *****green,
  mtype ****x,
  src_param *srcparam,
  det_param *detparam,
  phys_param *physparam
  /* Output */
);


int calc_phi(
     /*  Input */
  src_param *sources, /* Info about source locations and frequencies */
  det_param *dets,    /* Info about detector positions */
  mtype ****x,        /* mua and mus array: x[u,i,j,l] */
  phys_param *phys_param, /* Info about physical dimensions of problem */
    /* Output */
  mtype *****phi,     /* Calculated phi for each source: phi[k,i,j,l,c]*/
  mtype **fx         /* Calculated phi at measurement positions: fx[s,c]*/
 );

int calc_phi_fluo(
     /* Input */
  mtype ****fluo,    /* Fluorescent parameters, in complex form */
  mtype *****phi_x,  /* Photon fluence at excitation wavelength */
  src_param *sources, /* Info about source locations and frequencies */
  det_param *dets,    /* Info about detector positions */
  mtype ****x,        /* mua and mus array: x[u,i,j,l] */
  phys_param *phys_param, /* Info about physical dimensions of problem */
    /* Output */
  mtype *****phi_m,     /* Calculated phi for each source: phi[k,i,j,l,c]*/
  mtype **fx         /* Calculated phi at measurement positions: fx[s,c]*/
 );
 
int calc_green(
     /*  Input */
  src_param *sources, /* Info about source locations and frequencies */
  det_param *dets,    /* Info about detector positions */
  mtype ****x,        /* mua and mus array: x[u,i,j,l] */
  phys_param *phys_param, /* Info about physical dimensions of problem */
    /* Output */
  mtype *****green    /* Calculated green for each source: green[k,i,j,l,c]*/
 );
int calc_frechet_col_fluo_num(
  /* Input */
  int u, int i, int j, int l,
  mtype *****phi_x,
  mtype *****green,
  mtype ****fluo,
  src_param *sources,
  phys_param *phys_param,
    /* Output */
  mtype **At          /* At[s,c] */
) ;
int calc_frechet_col_fluo(
  /* Input */
  int u, int i, int j, int l,
  mtype *****phi_x,
  mtype *****green,
  mtype ****fluo,
  src_param *sources,
  phys_param *phys_param,
    /* Output */
  mtype **At          /* At[s,c] */
);
 
int calc_frechet_col(
  /* Input */
  int u, int i, int j, int l,
  mtype *****phi,
  mtype *****green,
  mtype ****x,
  src_param *sources,
  phys_param *phys_param,
    /* Output */
  mtype **At          /* At[s,c] */
);
int calc_frechet_row(
  /* Input */
  int s,
  mtype *****phi,
  mtype *****green,
  mtype ****x,
  src_param *srcparam,
  det_param *detparam,
  phys_param *physparam,
  /* Output */
  mtype *****At          /* At[u,i,j,l,c] */
);

int calc_frechet_row1(
  /* Input */
  int k,
  int m,
  mtype *****phi,
  mtype *****green,
  mtype ****x,
  src_param *srcparam,
  det_param *detparam,
  phys_param *physparam,
  /* Output */
  mtype ****At          /* At[u,i,j,l,c] */
);


/* Solves elliptic PDE in form of diffusion equation.  This is a
 * C function which uses parameters in a form similar to the rest
 * of the code. It calls fwd3df_ */
int fwdsolver_3d_fluo(
  /* Input */
  phys_param *phys_param, /* Info about physical dimensions of problem */
  src_param *sources, /* Info about source locations and frequencies */
  mtype ****fluo,    /* Fluorescent parameters, in complex form */
  mtype *****phi_x,  /* Photon fluence at excitation wavelength */
  int k,              /* Index describing which source to use in solution */
  mtype ****x,        /* mua and mus array: x[u,i,j,l] */
  int iguessflag,     /* =1 if initial guess is to be provided, =0 otherwise*/
  int ncycles,        /* Number of multigrid cycles to be used by MUDPACK */
  mtype ****phi       /* Input: initial guess, if iguessflag=1.  Should be
                       *   initialized to zeros otherwise.
                       * Output: solution of PDE.  In C, phi[i][j][l][c]
                       * In FORTRAN, phi(c+1,l+1,j+1,i+1) */
);  

int fwdsolver_3d(
  /* Input */
  phys_param *phys_param, /* Info about physical dimensions of problem */
  src_param *sources, /* Info about source locations and frequencies */
  int k,              /* Index describing which source to use in solution */
  mtype ****x,        /* mua and mus array: x[u,i,j,l] */
  int iguessflag,     /* =1 if initial guess is to be provided, =0 otherwise*/
  int ncycles,        /* Number of multigrid cycles to be used by MUDPACK */
  mtype ****phi       /* Input: initial guess, if iguessflag=1.  Should be
                       *   initialized to zeros otherwise.
                       * Output: solution of PDE.  In C, phi[i][j][l][c]
                       * In FORTRAN, phi(c+1,l+1,j+1,i+1) */
);  




/******************************************************************/
/***                        icd.c                               ***/
/******************************************************************/


prior_param *alloc_prior_param(
  double (*rho)(double, double, double, double), /* potential fn. rho(xi,xj,sigma,p) */
  dtype sigma, 
  dtype p, 
  int Nneighbor, 
  dtype *b);            /* coefficient for neighborhood relation */ 


void free_prior_param(prior_param *priorparam);


/* Calculate RMS error from yerror vector  */ 
double calc_yrmse(
  mtype ***yerror,
  phys_param *physparam);


double calc_rmse_with_lambda(
  mtype **yerror,
  mtype *lambda,
  phys_param *physparam);


/* Calculate lambda matrix */
void calc_lambda(
  mtype **y,             /* y[s][c]            */
  mtype **fx,            /* fx[s][c]           */ 
  phys_param *physparam,
  src_param *srcparam,
  det_param *detparam,
  mtype *lambda,    /* lambda[s] : output */
  mtype *alpha1);   /* output */ 


void ICD_params(
  mtype **At,          /* At[s,c]        */
  mtype **yerror,      /* yerror[s,c]    */
  mtype *lambda,       /* lambda[s]      */
  phys_param *physparam,
  src_param *srcparam,
  det_param *detparam,
  double *theta1,       /* output           */
  double *theta2);      /* output           */ 
/* Compute theta1 and theta2 for one-pixel ICD update                     */
/*  1) theta1 = -2 * Re{[col_{u,i,j,l}(At)]^H * rambda * error}           */
/*  2) theta2 = 2 * [col_{u,i,j,l}(At)]^H * rambda * [col_{u,i,j,l}(At)]  */	 
 

/* return the derivative of object funtion at x           */
double dev_obj_fn(
  mtype x,                       /* point being evaluated */ 
  mtype x_prev,                  /* value of last iteration */
  prior_param *priorparam,
  double theta1, 
  double theta2, 
  mtype *neighbor);


mtype obj_fn(
  phys_param *physparam,
  prior_param *prior_D,
  prior_param *prior_mua,
  mtype ****x,
  mtype alpha
) ;


/* find root of (derivative of object fn)=0 with a half-interval search           */
/* Solves equation (*f)(x,constants) = 0 on x in [a,b].                           */
/* In ODT, (*f)() will point dev_obj_fn().                                        */        
/* Requires that (*f)(a) and (*f)(b) have opposite signs.                         */
/* Returns code=0  if signs are opposite.                                         */
/* Returns code=1  if signs are both positive.                                    */
/* Returns code=-1 if signs are both negative.                                    */
double  root_find(
  double (*fp)(mtype, mtype, prior_param *, double, double, mtype *),
  double a,                /* minimum value of solution */
  double b,                /* maximum value of solution */
  double err,              /* accuarcy of solution */
  int *code,               /* error code */
  prior_param *priorparam, /* From this, the parameters are for (*f)() */
  mtype x_prev,            /* value of previous iteration */
  double theta1, 
  double theta2, 
  mtype *neighbor);        /* in 3D, neighbor[26] raster-ordered from [-1,-1,-1] to [1,1,1] */


/* store neighborhood pixel of x[u,i,j,l] into neighbor[] */
int get_neighbor(
  int u, int i, int j, int l,  /* center of the neighbor to be obtained */
  prior_param *priorparam,
  mtype ****x,                 /* x[u,i,j,l]   */
  mtype *neighbor);


/* return new value of x[u,i,j,l] obtained by one-pixel ICD  */
double ICD_pixel_update(
  int u, int i, int j, int l,  /* position of pixel being updated */
  mtype **At,                 /* At[c,s] for given u,i,j,l  */
  mtype ****x,                 /* x[u,i,j,l]   */
  mtype **yerror,             /* yerror[c,s]  */
  mtype *lambda,              /* lambda[s]   */
  prior_param *priorparam,
  phys_param *physparam,
  src_param *srcparam,
  det_param *detparam
);

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
)  ;



/* Update pixels performing one-cycle ICD iteration 
 * to find x^ = argmin_x { ||yerror-Ax||_rambda ^2 + x^t B x }.
 * One pixel update is performed by calling ICD_pixel_update().  
 * yerror is also updated. */
void ICD_update(
  mtype ****x,          /* initial value x[u][i][j][l] and output */
  mtype **yerror,      /* initial error and output               */   
  mtype *****phi,         /* phi[k][i][j][l][c]                     */ 
  mtype *****green,       /* green[m][i][j][l][c]                   */
  mtype **y,
  mtype **fx,
  mtype *lambda,        
  prior_param *priorparam_D, 
  prior_param *priorparam_mua, 
  src_param *srcparam,
  det_param *detparam,
  phys_param *physparam,
  config_param *config
);

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
  prior_param *priorparam_etamu,
  src_param *srcparam,
  det_param *detparam,
  phys_param *physparam,
  fluo_config_param *config
);


void rvssolver_sg_fluo(
  mtype ****fluo,              /* initial value of tau and etamuand output */
   mtype ****mu_m,
   mtype *****phi_x,
   mtype *****phi_m,
   mtype **meas_m,
  mtype **y,                /* observed phi, that is, desired value of f(x^) */
  phys_param *physparam,       /* physics parameters */
  src_param  *srcparam,        /* src_param array */
  det_param  *detparam,        /* detpos array    */
  prior_param *priorparam_tau,
  prior_param *priorparam_etamu,
  prior_param *priorparam_D,
  prior_param *priorparam_mua,
  mtype **yerror,             /* output : observ - f(x^)  */
  fluo_config_param *config);
 
          
/* Reverse solver using single-grid ICD algorithm.
 * This will be used in multigrid reverse solver subroutine rvssolver_mg().
 * Output : x (that is, mu_a & mu_s) and yerror */
void rvssolver_sg(
  mtype ****x,                 /* initial value of mu_a & mu_s and output */
  mtype **y,                   /* observed phi, that is, desired value of f(x^) */
  phys_param *physparam,       /* physics parameters */ 
  src_param  *srcparam,        /* src_param array */ 
  det_param  *detparam,        /* detpos array    */ 
  prior_param *priorparam_D,        
  prior_param *priorparam_mua,        
  mtype **yerror,
  config_param *config
);            /* output : observ - f(x^)  */


/* Reverse solver using multigrid will be considered in the future.
   This calls rvssolver_sg() */
void rvssolver_mg();


/******************************************************************/
/***                     icdsd.c                                ***/
/******************************************************************/


/* Calculate the derivatives of source weights w.r.t. mua and D.
   Used for updating source weights after each one pixel update.  */   
int calc_dsdx(
/* INPUT */
  int k,
  mtype **y,
  mtype **fx,
  mtype *lambda,
  mtype **At,                  /* At[s][c]  */
  phys_param *physparam,
  src_param *srcparam,
  det_param *detparam,
  config_param *config,
/* OUTPUT */
  mtype **dsdx                 /* dsdx[k][c] */   
);       


/* Calculate the derivatives of detector weights w.r.t. mua and D.
   Used for updating detector weights after each one pixel update.  */
int calc_dddx(
/* INPUT */
  int m,
  mtype **y,
  mtype **fx,
  mtype *lambda,
  mtype **At,                  /* At[s][c]  */
  phys_param *physparam,
  src_param *srcparam,
  det_param *detparam,
  config_param *config,
/* OUTPUT */
  mtype **dddx                 /* dddx[k][c] */
);


/* Update pixels performing one-cycle ICD iteration
   where the ds/dx and dd/dx are considered to calculate the Frechet derivatives dCost/dx.
   This routine does NOT perform the error correction for s-d calibration.
   (See ICD_update_sd2() . ) Actually this routine does not work well, 
   if a proper initialization for the s-d weights are not performed in rvs_solver_sd().
   x^ is found s.t.  x^ = argmin_x { ||yerror-Ax||_lambda ^2 + x^t B x }.
   One pixel update is performed by calling ICD_pixel_update().  */
void ICD_update_sd(
  mtype ****x,            /* initial value x[u][i][j][l] and output */
  mtype **yerror,         /* initial error and output               */
  mtype *****phi,         /* phi[k][i][j][l][c]                     */
  mtype *****green,       /* green[m][i][j][l][c]                   */
  mtype **y,
  mtype **fx,
  mtype *lambda,
  prior_param *priorparam_D,
  prior_param *priorparam_mua,
  src_param *srcparam,
  det_param *detparam,
  phys_param *physparam,
  config_param *config
);  


/* Update pixels performing one-cycle ICD iteration
   where the ds/dx and dd/dx are considered to calculate the Frechet derivatives dCost/dx.
   This routine perform the error correction for s-d calibration. (See ICD_update_sd(). )
   x^ is found s.t.  x^ = argmin_x { ||yerror-Ax||_lambda ^2 + x^t B x }.
   One pixel update is performed by calling ICD_pixel_update().  */
void ICD_update_sd2(
  mtype ****x,          /* initial value x[u][i][j][l] and output */
  mtype **yerror,      /* initial error and output               */
  mtype *****phi,         /* phi[k][i][j][l][c]                     */
  mtype *****green,       /* green[m][i][j][l][c]                   */
  mtype **y,
  mtype **fx,
  mtype *lambda,
  prior_param *priorparam_D,
  prior_param *priorparam_mua,
  src_param *srcparam,
  det_param *detparam,
  phys_param *physparam,
  config_param *config
);     

/******************************************************************/
/***                     homogeneous.c                          ***/
/******************************************************************/


/* Calculate the derivatives of global weights w.r.t. mua and D.
 * Used for calculate the gradient direction in homogeneous pixel update
 * with consideration of the derivatives of global weights w.r.t. mua and D. */
int calc_dwdx(
/* INPUT */
  mtype **y,
  mtype **fx,
  mtype *lambda,
  mtype ***SumAt,           /* SumAt[u][s][c] */
  phys_param *physparam,
  src_param *srcparam,
  det_param *detparam,
  config_param *config,
/* OUTPUT */
  double dwdx[2][2]);       /* dwdx[u][c] */   


/* Find the direction opposite to the gradient (dCost/dmua, dCost/dD),
   considering the effect of pixel update on the weights.
   This routine is originally designed to be used in Homo_update_Search2().   */
int calc_direction(
  mtype ****x,
  mtype **y,
  mtype **fx,
  mtype *****phi,
  mtype *****green,
  mtype *lambda,
  phys_param *physparam,
  src_param *srcparam,
  det_param *detparam,
  config_param *config,
  double dCdx[2]);          


/* Update pixels performing one-cycle Newton's method
   to find homogeneous x^ = argmin_x { ||yerror-Ax||_lambda ^2 + x^t B x }.
   Only one pixels  are assumed for each of mua and D. */
void Homo_update_Newton(
  mtype ****x,            /* initial value x[u][i][j][l] and output */
  mtype **yerror,         /* initial error and output               */
  mtype *****phi,         /* phi[k][i][j][l][c]                     */
  mtype *****green,       /* green[m][i][j][l][c]                   */
  mtype *lambda,
  src_param *srcparam,
  det_param *detparam,
  phys_param *physparam,
  config_param *config);


/* Update pixels performing one-cycle search method
   to find homogeneous x^ = argmin_x { ||yerror-Ax||_lambda ^2 + x^t B x }.
   Only one pixels  are assumed for each of mua and D. i
   In the directional search, the effect of pixel update on the weights
   are NOT considered, and search for mua and D is performed iteratively.
   (Compare this routine with Homo_update_Search2(). )  */
void Homo_update_Search(
  mtype ****x,
  mtype **y,
  mtype *lambda,
  mtype *****phi,         /* phi[k][i][j][l][c]                     */
  mtype *****green,       /* green[m][i][j][l][c]                   */
  mtype **yerror,
  src_param *srcparam,
  det_param *detparam,
  phys_param *physparam,
  config_param *config);         


/* Update pixels performing one-cycle search method
   to find homogeneous x^ = argmin_x { ||yerror-Ax||_lambda ^2 + x^t B x }.
   Only one pixels  are assumed for each of mua and D. i
   In the directional search, the effect of pixel update on the weights
   are considered, and search for mua and D is performed
   in the direction opposite to the gradient (dCost/dmua, dCost/dD).
   (Compare this routine with Homo_update_Search(). )  */
void Homo_update_Search2(
  mtype ****x,
  mtype **y,
  mtype **fx,
  mtype *lambda,
  mtype *****phi,         /* phi[k][i][j][l][c]                     */
  mtype *****green,       /* green[m][i][j][l][c]                   */
  mtype **yerror,
  src_param *srcparam,
  det_param *detparam,
  phys_param *physparam,
  config_param *config
);  



/******************************************************************/
/***                     calibrate.c                            ***/
/******************************************************************/

int findk(
  phys_param *physparam,
  int s);

int findm(
  phys_param *physparam,
  int s);

int one_src_calibrate(
  int k,
  mtype **y,
  mtype **fx,
  mtype *lambda,
  src_param *srcparam,
  det_param *detparam,
  phys_param *physparam);

int one_det_calibrate(
  int m,
  mtype **y,
  mtype **fx,
  mtype *lambda,
  src_param *srcparam,
  det_param *detparam,
  phys_param *physparam);

int src_det_calibrate(
  mtype **y,
  mtype **fx,
  mtype *lambda,
  src_param *srcparam,
  det_param *detparam,
  phys_param *physparam);

int global_weight_calibrate(
  mtype ****x,
  mtype **y,
  mtype **fx,
  mtype *lambda,
  src_param *srcparam,
  det_param *detparam,
  phys_param *physparam); 

int convert_fluo(
   int dir,
   phys_param *physparam,  /* physics parameters */
   src_param  *srcparam,   /* src_param array */
   mtype ****fluo_in,
   mtype ****fluo_out
);



int convert_config(fluo_config_param *config, config_param *config_x);


mtype obj_fn_fluo(
  phys_param *physparam,
  prior_param *prior_tau,
  prior_param *prior_etamu,
  prior_param *prior_D,
  prior_param *prior_mua,
  mtype ****fluo_lin,
  mtype ****mu_m,
  mtype alpha
);
int update_frechet_col_fluo(
   int i, int j, int l, int u,
   mtype *****phi_x,
   mtype *****green,
   mtype ****fluo, 
   src_param *sources, 
   phys_param *physparam, 
   mtype **At);

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
);


int findr(src_param *src, int karg);


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
);


mtype  bracket(
  mtype as,
  mtype bs,
  int N,
  arglist *args,
  mtype (*f)(arglist *, mtype )
);

mtype golden(mtype ax,mtype bx,mtype cx,mtype (*f)(arglist *, mtype),
  arglist *args,
  mtype tol,int itmax, mtype *xmin);


mtype ff( arglist *args , mtype x);

#define CALI_FLAG 0
#if CALI_FLAG
  #define calimultr(a_r, a_i, srcparam, detparam, physparam,k,m) {    \
        cplxmultr( \
	   cplxmultr(a_r, a_i, \
     	     cplxmultr(srcparam[k].calir, srcparam[k].calii, detparam[m].calir, detparam[m].calii), \
	     cplxmulti(srcparam[k].calir, srcparam[k].calii, detparam[m].calir, detparam[m].calii)  \
           ), \
	   cplxmulti(a_r, a_i,  \
             cplxmultr(srcparam[k].calir, srcparam[k].calii, detparam[m].calir, detparam[m].calii), \
	     cplxmulti(srcparam[k].calir, srcparam[k].calii, detparam[m].calir, detparam[m].calii)  \
           ), \
           physparam->wgtr, physparam->wgti \
	 )
#else
  #define calimultr(a_r, a_i, srcparam, detparam, physparam,k,m) a_r   
#endif

#if CALI_FLAG
  #define calimulti(a_r, a_i, srcparam, detparam, physparam,k,m) {    \
        cplxmulti( \
	   cplxmultr(a_r, a_i, \
     	     cplxmultr(srcparam[k].calir, srcparam[k].calii, detparam[m].calir, detparam[m].calii), \
	     cplxmulti(srcparam[k].calir, srcparam[k].calii, detparam[m].calir, detparam[m].calii)\
           ), \
	   cplxmulti(a_r, a_i, \
             cplxmultr(srcparam[k].calir, srcparam[k].calii, detparam[m].calir, detparam[m].calii), \
	     cplxmulti(srcparam[k].calir, srcparam[k].calii, detparam[m].calir, detparam[m].calii)  \
           ), \
	   physparam->wgtr, physparam->wgti \
	 )
#else
  #define calimulti(a_r, a_i, srcparam, detparam, physparam,k,m) a_i   
#endif

