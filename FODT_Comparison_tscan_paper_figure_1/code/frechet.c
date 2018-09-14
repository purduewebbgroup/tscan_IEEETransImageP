/*
 * Compare Frechet derivative computed by Jong's algorithm
 * with the measured derivative 
 * 
 * HISTORY: 
 *
 * Seungseok Oh 
 * 10/7/2000
 *
 * Adam Milstein
 * 10/11/2000
 * Fixed sign error in Green's function.
 * Added 3 arguments to fwdsolverf_ to allow workspace
 *   size, a flag for indicating initial guess or not,
 *   and the number of desired multigrid cycles.
 *
 * 10/14/2000
 * In 3 places, I fixed it so that i goes from 1 to nnx and
 *  j goes from 1 to nny, instead of something else, which 
 *  would be wrong.  This slipped by us before because nnx 
 *  was equal to nny for our test cases!
 * 
 */


#include "defs.h"

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
  mtype ****phi_m     /* Input: initial guess, if iguessflag=1.  Should be
                       *   initialized to zeros otherwise.
                       * Output: solution of PDE.  In C, phi[i][j][l][c]
                       * In FORTRAN, phi(c+1,l+1,j+1,i+1) */
)
{
  int i,j,l, mudpkworksize;
  int nnx, nny, nnz;
  dtype dx,dy,dz;
  mtype xmin, xmax, ymin, ymax,zmin, zmax;
  mtype fluo_i, fluo_r;
  mtype tau, eta_mu;
  FILE *fp=NULL; /* DEBUG */
  mtype ***rhsworkspace, ***phiworkspace;
  mtype *mudpkworkspace;

  dtype omega,v;
  mtype ****alpha, ****beta,  
        ****rhs;

  nnx=phys_param->Ni;
  nny=phys_param->Nj;
  nnz=phys_param->Nl;
  
  xmin=(mtype)(phys_param->xmin);
  xmax=(mtype)(phys_param->xmax);
  ymin=(mtype)(phys_param->ymin);
  ymax=(mtype)(phys_param->ymax);
  zmin=(mtype)(phys_param->zmin);
  zmax=(mtype)(phys_param->zmax);
  
  dx=(xmax-xmin)/((dtype)(nnx-1));
  dy=(ymax-ymin)/((dtype)(nny-1));
  dz=(zmax-zmin)/((dtype)(nnz-1));

  /* allocate and initialize matrices for PDE coefficients  */

  alpha= (mtype ****)multialloc(sizeof(mtype), 4, nnx, nny, nnz, 2);
  beta = (mtype ****)multialloc(sizeof(mtype), 4, nnx, nny, nnz, 2);
  rhs  = (mtype ****)multialloc(sizeof(mtype), 4, nnx, nny, nnz, 2);

/* Stuff for interfacing with Mudpack */

  rhsworkspace  = (mtype ***)multialloc(2*sizeof(mtype), 3, nnz, nny, nnx);
  phiworkspace  = (mtype ***)multialloc(2*sizeof(mtype), 3, nnz, nny, nnx);


/*  This number changes if you are using a MUDPACK solver other than
    cud3.f .  See the MUDPACK documentation (cud3.d, in particular) for
    more details */

  mudpkworksize = 3*(nnx+2)*(nny+2)*(nnz+2)*(10+0+0+0);
  mudpkworkspace  = (mtype *)malloc(mudpkworksize*sizeof(mtype));
  omega=sources[k].omega;
  v = phys_param-> v;
 
  for(i=0; i<nnx; i++)
  for(j=0; j<nny; j++)
  for(l=0; l<nnz; l++)
  {
      tau=fluo[0][i][j][l];
      eta_mu=fluo[1][i][j][l];
      fluo_r=eta_mu/(1.0+omega*omega*tau*tau); 
      fluo_i=-eta_mu*omega*tau/(1.0+omega*omega*tau*tau);
      rhs[i][j][l][0]=-(mtype)cplxmultr(phi_x[k][i][j][l][0], phi_x[k][i][j][l][1],
                                        fluo_r, fluo_i);
      rhs[i][j][l][1]=-(mtype)cplxmulti(phi_x[k][i][j][l][0], phi_x[k][i][j][l][1],
                                        fluo_r, fluo_i);
      if (iguessflag==0) 
      {
         phi_m[i][j][l][0]=(mtype)(0.0);
         phi_m[i][j][l][1]=(mtype)(0.0);
      } 
  }


  

  for(i=0; i<nnx; i++)
  for(j=0; j<nny; j++)
  for(l=0; l<nnz; l++)
  {
      alpha[i][j][l][0] = x[1][i][j][l]; 
      alpha[i][j][l][1] = 0.0; 
      beta[i][j][l][0]  = -x[0][i][j][l];
      beta[i][j][l][1]  = -omega/v;
  }


  fwd3df_(&nnx,&nny,&nnz,
          &xmin,&xmax,&ymin,&ymax,&zmin,&zmax,
          &alpha[0][0][0][0], 
          &beta[0][0][0][0], 
          &rhs[0][0][0][0], 
          &iguessflag, 
          &mudpkworkspace[0],
          &mudpkworksize,
          &rhsworkspace[0][0][0], 
          &phiworkspace[0][0][0], 
          &ncycles,  
          &phi_m[0][0][0][0]);
/*
     fp = datOpen("rhs.dat", "w+b");
     write_float_array(fp, "rhs", &rhs[0][0][0][0], 4, nnx, nny, nnz, 2);
     datClose(fp);

     fp = datOpen("phi.dat", "w+b");
     write_float_array(fp, "phi", &phi_m[0][0][0][0], 4, nnx, nny, nnz, 2);
     datClose(fp);
*/
  multifree(alpha,4);
  multifree(beta,4);
  multifree(rhs,4);
  multifree(rhsworkspace,3);
  multifree(phiworkspace,3);
  free(mudpkworkspace);

  return(0); 
} 



/* Solves elliptic PDE in form of diffusion equation.  This is a
 * C function which uses parameters in a form similar to the rest
 * of the code. It calls fwd3df_ */
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
)
{
  int i,j,l, mudpkworksize;
  int i_l, i_h, j_l, j_h, l_l, l_h;
  int nnx, nny, nnz;
  dtype dx,dy,dz;
  mtype xmin, xmax, ymin, ymax,zmin, zmax;
  double x_l, x_h, y_l, y_h, z_l, z_h;
  double x_s, y_s, z_s;
  double src_amplitude;

  mtype ***rhsworkspace, ***phiworkspace;
  mtype *mudpkworkspace;

  dtype omega,v;
  mtype ****alpha, ****beta,  
        ****rhs;

  nnx=phys_param->Ni;
  nny=phys_param->Nj;
  nnz=phys_param->Nl;
  
  xmin=(mtype)(phys_param->xmin);
  xmax=(mtype)(phys_param->xmax);
  ymin=(mtype)(phys_param->ymin);
  ymax=(mtype)(phys_param->ymax);
  zmin=(mtype)(phys_param->zmin);
  zmax=(mtype)(phys_param->zmax);
  
  dx=(xmax-xmin)/((dtype)(nnx-1));
  dy=(ymax-ymin)/((dtype)(nny-1));
  dz=(zmax-zmin)/((dtype)(nnz-1));

  /* allocate and initialize matrices for PDE coefficients  */

  alpha= (mtype ****)multialloc(sizeof(mtype), 4, nnx, nny, nnz, 2);
  beta = (mtype ****)multialloc(sizeof(mtype), 4, nnx, nny, nnz, 2);
  rhs  = (mtype ****)multialloc(sizeof(mtype), 4, nnx, nny, nnz, 2);

/* Stuff for interfacing with Mudpack */

  rhsworkspace  = (mtype ***)multialloc(2*sizeof(mtype), 3, nnz, nny, nnx);
  phiworkspace  = (mtype ***)multialloc(2*sizeof(mtype), 3, nnz, nny, nnx);


/*  This number changes if you are using a MUDPACK solver other than
    cud3.f .  See the MUDPACK documentation (cud3.d, in particular) for
    more details */

  mudpkworksize = 3*(nnx+2)*(nny+2)*(nnz+2)*(10+0+0+0);
  mudpkworkspace  = (mtype *)malloc(mudpkworksize*2*sizeof(mtype));
 
  for(i=0; i<nnx; i++)
  for(j=0; j<nny; j++)
  for(l=0; l<nnz; l++)
  {
      rhs[i][j][l][0]=(mtype)(0.0);
      rhs[i][j][l][1]=(mtype)(0.0); 
      if (iguessflag==0) 
      {
         phi[i][j][l][0]=(mtype)(0.0);
         phi[i][j][l][1]=(mtype)(0.0);
      } 
  }

  x_s=sources[k].x;
  y_s=sources[k].y;
  z_s=sources[k].z;
  omega=sources[k].omega;
  v = phys_param-> v;

  i_l=(int)floor((x_s-xmin)/dx);
  j_l=(int)floor((y_s-ymin)/dy);
  l_l=(int)floor((z_s-zmin)/dz);
  x_l=xmin+dx*(double)(i_l);
  y_l=ymin+dy*(double)(j_l);
  z_l=zmin+dz*(double)(l_l);
  i_h=i_l+1;
  j_h=j_l+1;
  l_h=l_l+1;
  x_h=x_l+dx;
  y_h=y_l+dy;
  z_h=z_l+dz;

  src_amplitude= (-1.0*sources[k].beta/(dx*dy*dz));
  rhs[i_h][j_h][l_h][0] = (mtype)(src_amplitude*(x_s-x_l)/dx*(y_s-y_l)/dy*(z_s-z_l)/dz);
  rhs[i_h][j_h][l_l][0] = (mtype)(src_amplitude*(x_s-x_l)/dx*(y_s-y_l)/dy*(z_h-z_s)/dz);
  rhs[i_h][j_l][l_h][0] = (mtype)(src_amplitude*(x_s-x_l)/dx*(y_h-y_s)/dy*(z_s-z_l)/dz);
  rhs[i_h][j_l][l_l][0] = (mtype)(src_amplitude*(x_s-x_l)/dx*(y_h-y_s)/dy*(z_h-z_s)/dz);
  
  rhs[i_l][j_h][l_h][0] = (mtype)(src_amplitude*(x_h-x_s)/dx*(y_s-y_l)/dy*(z_s-z_l)/dz);
  rhs[i_l][j_h][l_l][0] = (mtype)(src_amplitude*(x_h-x_s)/dx*(y_s-y_l)/dy*(z_h-z_s)/dz);
  rhs[i_l][j_l][l_h][0] = (mtype)(src_amplitude*(x_h-x_s)/dx*(y_h-y_s)/dy*(z_s-z_l)/dz);
  rhs[i_l][j_l][l_l][0] = (mtype)(src_amplitude*(x_h-x_s)/dx*(y_h-y_s)/dy*(z_h-z_s)/dz);


  for(i=0; i<nnx; i++)
  for(j=0; j<nny; j++)
  for(l=0; l<nnz; l++)
  {
      alpha[i][j][l][0] = x[1][i][j][l]; 
      alpha[i][j][l][1] = 0.0; 
      beta[i][j][l][0]  = -x[0][i][j][l];
      beta[i][j][l][1]  = -omega/v;
  }

  fwd3df_(&nnx,&nny,&nnz,
          &xmin,&xmax,&ymin,&ymax,&zmin,&zmax,
          &alpha[0][0][0][0], 
          &beta[0][0][0][0], 
          &rhs[0][0][0][0], 
          &iguessflag, 
          &mudpkworkspace[0],
          &mudpkworksize,
          &rhsworkspace[0][0][0], 
          &phiworkspace[0][0][0], 
          &ncycles,  
          &phi[0][0][0][0]);

  multifree(alpha,4);
  multifree(beta,4);
  multifree(rhs,4);
  multifree(rhsworkspace,3);
  multifree(phiworkspace,3);
  free(mudpkworkspace);

  return(0); 
} 

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
 )
{
  int k,m,s, i_l,j_l,l_l, i_h, j_h, l_h ;
  int nnx, nny, nnz;
  LINK2D sparse, src_element;
  LINK   det_list;
  dtype  xmin, xmax, ymin, ymax, zmin, zmax, dx, dy, dz;
  dtype  x_l, x_h, y_l, y_h, z_l, z_h, x_d, y_d, z_d;

  sparse=phys_param->sparse;
  
  s=0;
  for(k=0; k<phys_param->K; k++){

    fwdsolver_3d_fluo( phys_param, sources, fluo, phi_x, 
       k, x, 0, DEFAULT_CYCLES, phi_m[k]);


    src_element=find_element_2d(sparse, k);
    if(src_element==NULL){
      fprintf(stderr, "calc_phi: problem with linked list. Aborting.\n");
      return -1;
    }
    det_list=src_element->d;

    nnx=phys_param->Ni;
    nny=phys_param->Nj;
    nnz=phys_param->Nl;
    xmin=(phys_param->xmin);
    xmax=(phys_param->xmax);
    ymin=(phys_param->ymin);
    ymax=(phys_param->ymax);
    zmin=(phys_param->zmin);
    zmax=(phys_param->zmax);
    
    dx=(xmax-xmin)/((dtype)(nnx-1));
    dy=(ymax-ymin)/((dtype)(nny-1));
    dz=(zmax-zmin)/((dtype)(nnz-1));

    for(m=0; m<phys_param->M; m++){
      x_d=dets[m].x;
      y_d=dets[m].y;
      z_d=dets[m].z;

      i_l=(int)floor((x_d-xmin)/dx);
      j_l=(int)floor((y_d-ymin)/dy);
      l_l=(int)floor((z_d-zmin)/dz);
      x_l=xmin+dx*(double)(i_l);
      y_l=ymin+dy*(double)(j_l);
      z_l=zmin+dz*(double)(l_l);
      i_h=i_l+1;
      j_h=j_l+1;
      l_h=l_l+1;
      x_h=x_l+dx;
      y_h=y_l+dy;
      z_h=z_l+dz;


      if(find_data_in_list(det_list,m) >= 0){
        fx[s][0]=(x_d-x_l)/dx*(y_d-y_l)/dy*(z_d-z_l)/dz*phi_m[k][i_h][j_h][l_h][0]
                +(x_d-x_l)/dx*(y_d-y_l)/dy*(z_h-z_d)/dz*phi_m[k][i_h][j_h][l_l][0]
                +(x_d-x_l)/dx*(y_h-y_d)/dy*(z_d-z_l)/dz*phi_m[k][i_h][j_l][l_h][0]
                +(x_d-x_l)/dx*(y_h-y_d)/dy*(z_h-z_d)/dz*phi_m[k][i_h][j_l][l_l][0]
                +(x_h-x_d)/dx*(y_d-y_l)/dy*(z_d-z_l)/dz*phi_m[k][i_l][j_h][l_h][0]
                +(x_h-x_d)/dx*(y_d-y_l)/dy*(z_h-z_d)/dz*phi_m[k][i_l][j_h][l_l][0]
                +(x_h-x_d)/dx*(y_h-y_d)/dy*(z_d-z_l)/dz*phi_m[k][i_l][j_l][l_h][0]
                +(x_h-x_d)/dx*(y_h-y_d)/dy*(z_h-z_d)/dz*phi_m[k][i_l][j_l][l_l][0];
 
        fx[s][1]=(x_d-x_l)/dx*(y_d-y_l)/dy*(z_d-z_l)/dz*phi_m[k][i_h][j_h][l_h][1]
                +(x_d-x_l)/dx*(y_d-y_l)/dy*(z_h-z_d)/dz*phi_m[k][i_h][j_h][l_l][1]
                +(x_d-x_l)/dx*(y_h-y_d)/dy*(z_d-z_l)/dz*phi_m[k][i_h][j_l][l_h][1]
                +(x_d-x_l)/dx*(y_h-y_d)/dy*(z_h-z_d)/dz*phi_m[k][i_h][j_l][l_l][1]
                +(x_h-x_d)/dx*(y_d-y_l)/dy*(z_d-z_l)/dz*phi_m[k][i_l][j_h][l_h][1]
                +(x_h-x_d)/dx*(y_d-y_l)/dy*(z_h-z_d)/dz*phi_m[k][i_l][j_h][l_l][1]
                +(x_h-x_d)/dx*(y_h-y_d)/dy*(z_d-z_l)/dz*phi_m[k][i_l][j_l][l_h][1]
                +(x_h-x_d)/dx*(y_h-y_d)/dy*(z_h-z_d)/dz*phi_m[k][i_l][j_l][l_l][1];
                
        s++;
      }

    }
  }
  return 0;
}

int calc_phi(
     /*  Input */
  src_param *sources, /* Info about source locations and frequencies */
  det_param *dets,    /* Info about detector positions */
  mtype ****x,        /* mua and mus array: x[u,i,j,l] */
  phys_param *phys_param, /* Info about physical dimensions of problem */
    /* Output */
  mtype *****phi,     /* Calculated phi for each source: phi[k,i,j,l,c]*/
  mtype **fx         /* Calculated phi at measurement positions: fx[s,c]*/
 )
{
  int k,m,s, i_l,j_l,l_l, i_h, j_h, l_h ;
  int nnx, nny, nnz;
  LINK2D sparse, src_element;
  LINK   det_list;
  dtype  xmin, xmax, ymin, ymax, zmin, zmax, dx, dy, dz;
  dtype  x_l, x_h, y_l, y_h, z_l, z_h, x_d, y_d, z_d;

  sparse=phys_param->sparse;
  /****************/
  printf("calc_phi\n");
  /*******************/
  s=0;
  for(k=0; k<phys_param->K; k++){
    fwdsolver_3d(phys_param, sources,  
       k,x,0, DEFAULT_CYCLES, phi[k]);

    src_element=find_element_2d(sparse, k);
    if(src_element==NULL){
      fprintf(stderr, "calc_phi: problem with linked list. Aborting.\n");
      return -1;
    }
    det_list=src_element->d;

    nnx=phys_param->Ni;
    nny=phys_param->Nj;
    nnz=phys_param->Nl;
    xmin=(phys_param->xmin);
    xmax=(phys_param->xmax);
    ymin=(phys_param->ymin);
    ymax=(phys_param->ymax);
    zmin=(phys_param->zmin);
    zmax=(phys_param->zmax);
    
    dx=(xmax-xmin)/((dtype)(nnx-1));
    dy=(ymax-ymin)/((dtype)(nny-1));
    dz=(zmax-zmin)/((dtype)(nnz-1));

    for(m=0; m<phys_param->M; m++){
      x_d=dets[m].x;
      y_d=dets[m].y;
      z_d=dets[m].z;

      i_l=(int)floor((x_d-xmin)/dx);
      j_l=(int)floor((y_d-ymin)/dy);
      l_l=(int)floor((z_d-zmin)/dz);
      x_l=xmin+dx*(double)(i_l);
      y_l=ymin+dy*(double)(j_l);
      z_l=zmin+dz*(double)(l_l);
      i_h=i_l+1;
      j_h=j_l+1;
      l_h=l_l+1;
      x_h=x_l+dx;
      y_h=y_l+dy;
      z_h=z_l+dz;


      if(find_data_in_list(det_list,m) >= 0){
        fx[s][0]=(x_d-x_l)/dx*(y_d-y_l)/dy*(z_d-z_l)/dz*phi[k][i_h][j_h][l_h][0]
                +(x_d-x_l)/dx*(y_d-y_l)/dy*(z_h-z_d)/dz*phi[k][i_h][j_h][l_l][0]
                +(x_d-x_l)/dx*(y_h-y_d)/dy*(z_d-z_l)/dz*phi[k][i_h][j_l][l_h][0]
                +(x_d-x_l)/dx*(y_h-y_d)/dy*(z_h-z_d)/dz*phi[k][i_h][j_l][l_l][0]
                +(x_h-x_d)/dx*(y_d-y_l)/dy*(z_d-z_l)/dz*phi[k][i_l][j_h][l_h][0]
                +(x_h-x_d)/dx*(y_d-y_l)/dy*(z_h-z_d)/dz*phi[k][i_l][j_h][l_l][0]
                +(x_h-x_d)/dx*(y_h-y_d)/dy*(z_d-z_l)/dz*phi[k][i_l][j_l][l_h][0]
                +(x_h-x_d)/dx*(y_h-y_d)/dy*(z_h-z_d)/dz*phi[k][i_l][j_l][l_l][0];
 
        fx[s][1]=(x_d-x_l)/dx*(y_d-y_l)/dy*(z_d-z_l)/dz*phi[k][i_h][j_h][l_h][1]
                +(x_d-x_l)/dx*(y_d-y_l)/dy*(z_h-z_d)/dz*phi[k][i_h][j_h][l_l][1]
                +(x_d-x_l)/dx*(y_h-y_d)/dy*(z_d-z_l)/dz*phi[k][i_h][j_l][l_h][1]
                +(x_d-x_l)/dx*(y_h-y_d)/dy*(z_h-z_d)/dz*phi[k][i_h][j_l][l_l][1]
                +(x_h-x_d)/dx*(y_d-y_l)/dy*(z_d-z_l)/dz*phi[k][i_l][j_h][l_h][1]
                +(x_h-x_d)/dx*(y_d-y_l)/dy*(z_h-z_d)/dz*phi[k][i_l][j_h][l_l][1]
                +(x_h-x_d)/dx*(y_h-y_d)/dy*(z_d-z_l)/dz*phi[k][i_l][j_l][l_h][1]
                +(x_h-x_d)/dx*(y_h-y_d)/dy*(z_h-z_d)/dz*phi[k][i_l][j_l][l_l][1];
                
        s++;
      }

    }
  }
  return 0;
}

int add_detector_noise(phys_param *physparam, dtype alpha_fixed, 
                                              dtype max_snr, /* Milstein */
                                              mtype **meas, mtype *snr){
  int s;
  double noisymeas_real, noisymeas_imag, mag;
  double alpha_s;
  srandom2(time(NULL));

  for(s=0; s<physparam->S; s++){
    mag=sqrt(AbsSquare(meas[s][0],meas[s][1]));

    if(alpha_fixed>1.0e-15){
      snr[s]=1./alpha_fixed*mag;
    }
    else{
      snr[s]=1.0e15;
    }
    alpha_s=alpha_fixed;
    if(max_snr<100.0 &&  10.0*log10(snr[s]) > max_snr ){
      alpha_s=mag/pow(10.0, max_snr/10.0);
      snr[s]=1./alpha_s*mag;
    }

    noisymeas_real=meas[s][0]+sqrt(.5*alpha_s*mag)*normal();
    noisymeas_imag=meas[s][1]+sqrt(.5*alpha_s*mag)*normal();

    meas[s][0]=noisymeas_real;
    meas[s][1]=noisymeas_imag;
  }
  return 0; 
}




/******
int calc_gradient(
  mtype ****field,
  mtype **grad,
  int i, int j, int l,
  dtype dx, dtype dy, dtype dz,
  int opt)
{
  int c, x, y;
  dtype h[2][3][3]
    = { { {0, 0, 0}, {0,0.5,0},{0,0,0} },
        { {1./18., 1./18., 1./18}, {1./18., 1./18., 1./18}, {1./18., 1./18., 1./18} } };

  for (c=0; c<2; c++) {
    grad[0][c] = 0.0;
    grad[1][c] = 0.0;
    grad[2][c] = 0.0;

    for (x=0; x<3; x++)
    for (y=0; y<3; y++) {
      grad[0][c] += (-field[i-1][j+x-1][l+y-1][c] + field[i+1][j+x-1][l+y-1][c]) * h[opt][x][y] / (dx * 2.);
      grad[1][c] += (-field[i+x-1][j-1][l+y-1][c] + field[i+x-1][j+1][l+y-1][c]) * h[opt][x][y] / (dy * 2.);
      grad[2][c] += (-field[i+x-1][j+y-1][l-1][c] + field[i+x-1][j+y-1][l+1][c]) * h[opt][x][y] / (dz * 2.);
    }
  }

  return 0;
}
****/

int calc_gradient(
  mtype ****field,
  mtype **grad,
  int i, int j, int l,
  dtype dx, dtype dy, dtype dz,
  int opt)
{
  int c;

  for (c=0; c<2; c++) {
    grad[0][c] = (-field[i-1][j][l][c] + field[i+1][j][l][c]) / (dx*2.);
    grad[1][c] = (-field[i][j-1][l][c] + field[i][j+1][l][c]) / (dy*2.);
    grad[2][c] = (-field[i][j][l-1][c] + field[i][j][l+1][c]) / (dz*2.);
  }

  return 0;
}

int calc_green(
     /*  Input */
  src_param *sources, /* Info about source locations and frequencies */
  det_param *dets,    /* Info about detector positions */
  mtype ****x,        /* mua and mus array: x[u,i,j,l] */
  phys_param *phys_param, /* Info about physical dimensions of problem */
    /* Output */
  mtype *****green    /* Calculated green for each source: green[k,i,j,l,c]*/
 )
{
  int m;
  src_param *green_src_at_det_pos=NULL;
  double x_d, y_d, z_d;

  green_src_at_det_pos=(src_param *)malloc(1*sizeof(src_param));
  printf("calc_green\n");

  for(m=0; m<phys_param->M; m++){
    x_d = dets[m].x;
    y_d = dets[m].y;
    z_d = dets[m].z;

    green_src_at_det_pos->x=x_d;
    green_src_at_det_pos->y=y_d;
    green_src_at_det_pos->z=z_d;
    green_src_at_det_pos->omega=dets[m].omega;
    green_src_at_det_pos->beta=(dtype)1.0;
  
    fwdsolver_3d(phys_param, green_src_at_det_pos, 
       0,x,0, DEFAULT_CYCLES, green[m]);

  }
  free(green_src_at_det_pos);
  return 0;
}

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
)
{
  int k,m,s;
  int nnx, nny, nnz;
  dtype xmin, xmax, ymin, ymax, zmin, zmax, dx, dy, dz, A;
  dtype D, mua;
  dtype v, omega;
  dtype greenr, greeni, phir, phii, tempr, tempi;
  mtype **gradgreen, **gradphi;    /* gradient[x/y/z][r/i] */
  double dgr, dgi, dpr, dpi;

  LINK2D sparse, src_element;
  LINK   det_list;

  nnx=phys_param->Ni;
  nny=phys_param->Nj;
  nnz=phys_param->Nl;

  v=phys_param->v;
  xmin=phys_param->xmin;
  xmax=phys_param->xmax;
  ymin=phys_param->ymin;
  ymax=phys_param->ymax;
  zmin=phys_param->zmin;
  zmax=phys_param->zmax;
 
  dx=(xmax-xmin)/((dtype)(nnx-1));
  dy=(ymax-ymin)/((dtype)(nny-1));
  dz=(zmax-zmin)/((dtype)(nnz-1));
  A = dx * dy * dz;

  gradgreen = multialloc(sizeof(mtype), 2, 3, 2);
  gradphi   = multialloc(sizeof(mtype), 2, 3, 2);

  sparse=phys_param->sparse;
  s=0;
  for(k=0; k<phys_param->K; k++){
    src_element=find_element_2d(sparse, k);
    if(src_element==NULL){
      fprintf(stderr, "calc_frechet_col: problem with linked list. Aborting.\n");
      return -1;
    }
    det_list=src_element->d;

    for(m=0; m<phys_param->M; m++) {
      if(find_data_in_list(det_list,m) >= 0){
        omega=sources[k].omega;
        mua=(dtype)x[0][i][j][l];
        D  =(dtype)x[1][i][j][l];
        greenr=green[m][i][j][l][0];
        greeni=green[m][i][j][l][1];

        phir=phi[k][i][j][l][0];
        phii=phi[k][i][j][l][1];

        tempr = (greenr * phir - greeni * phii) * A;
        tempi = (greenr * phii + greeni * phir) * A;

        if(u==MUA) {
          At[s][0] = -tempr;
          At[s][1] = -tempi;
        }

        else if(u==MUS) {
/*  In case of finite-element-like gradient,
          calc_gradient(green[m],gradgreen,i,j,l,dx,dy,dz,0);
          calc_gradient(phi[k],  gradphi,  i,j,l,dx,dy,dz,0);

          At[s][0] = 0.0;
          At[s][1] = 0.0;

          for (ii=0; ii<3; ii++) {
            At[s][0] += gradgreen[ii][0] * gradphi[ii][0]
                         - gradgreen[ii][1] * gradphi[ii][1];
            At[s][1] += gradgreen[ii][0] * gradphi[ii][1]
                         + gradgreen[ii][1] * gradphi[ii][0];
          }

          At[s][0] *= -A; 
          At[s][1] *= -A;
*/
/* In case of Jong's Formular-like method ,
          cDr = mua/D;
          cDi = omega/v/D;
          At[s][0] = tempr * cDr - tempi * cDi;
          At[s][1] = tempr * cDi + tempi * cDr;
*/

/* In case of Professors' method, */
          At[s][0] = 0.0;
          At[s][1] = 0.0;

          dgr = (green[m][i][j][l][0] - green[m][i-1][j][l][0])/dx;
          dgi = (green[m][i][j][l][1] - green[m][i-1][j][l][1])/dx;
          dpr = (phi[k][i][j][l][0] - phi[k][i-1][j][l][0])/dx;
          dpi = (phi[k][i][j][l][1] - phi[k][i-1][j][l][1])/dx;
          At[s][0] += dgr * dpr - dgi * dpi;
          At[s][1] += dgi * dpr + dgr * dpi;

          dgr = (green[m][i+1][j][l][0] - green[m][i][j][l][0])/dx;
          dgi = (green[m][i+1][j][l][1] - green[m][i][j][l][1])/dx;
          dpr = (phi[k][i+1][j][l][0] - phi[k][i][j][l][0])/dx;
          dpi = (phi[k][i+1][j][l][1] - phi[k][i][j][l][1])/dx;
          At[s][0] += dgr * dpr - dgi * dpi;
          At[s][1] += dgi * dpr + dgr * dpi;

          dgr = (green[m][i][j][l][0] - green[m][i][j-1][l][0])/dy;
          dgi = (green[m][i][j][l][1] - green[m][i][j-1][l][1])/dy;
          dpr = (phi[k][i][j][l][0] - phi[k][i][j-1][l][0])/dy;
          dpi = (phi[k][i][j][l][1] - phi[k][i][j-1][l][1])/dy;
          At[s][0] += dgr * dpr - dgi * dpi;
          At[s][1] += dgi * dpr + dgr * dpi;

          dgr = (green[m][i][j+1][l][0] - green[m][i][j][l][0])/dy;
          dgi = (green[m][i][j+1][l][1] - green[m][i][j][l][1])/dy;
          dpr = (phi[k][i][j+1][l][0] - phi[k][i][j][l][0])/dy;
          dpi = (phi[k][i][j+1][l][1] - phi[k][i][j][l][1])/dy;
          At[s][0] += dgr * dpr - dgi * dpi;
          At[s][1] += dgi * dpr + dgr * dpi;

          dgr = (green[m][i][j][l][0] - green[m][i][j][l-1][0])/dz;
          dgi = (green[m][i][j][l][1] - green[m][i][j][l-1][1])/dz;
          dpr = (phi[k][i][j][l][0] - phi[k][i][j][l-1][0])/dz;
          dpi = (phi[k][i][j][l][1] - phi[k][i][j][l-1][1])/dz;
          At[s][0] += dgr * dpr - dgi * dpi;
          At[s][1] += dgi * dpr + dgr * dpi;

          dgr = (green[m][i][j][l+1][0] - green[m][i][j][l][0])/dz;
          dgi = (green[m][i][j][l+1][1] - green[m][i][j][l][1])/dz;
          dpr = (phi[k][i][j][l+1][0] - phi[k][i][j][l][0])/dz;
          dpi = (phi[k][i][j][l+1][1] - phi[k][i][j][l][1])/dz;
        
          At[s][0] += dgr * dpr - dgi * dpi;
          At[s][1] += dgi * dpr + dgr * dpi;

          At[s][0] *= -A / 2.0;
          At[s][1] *= -A / 2.0;
/**/
        }
        else {
/*          fprintf(stderr, "ERROR! Invalid u!\n"); */
           printf("ERROR! Invalid u!\n"); 
          exit(1);
        }
        s++;
      }
    }/* End for m */
  }/* End for k*/
  multifree(gradgreen, 2);
  multifree(gradphi, 2);

  return 0;
}
int update_frechet_col_fluo(
   int i, int j, int l, int u,
   mtype *****phi_x,
   mtype *****green,
   mtype ****fluo,
   src_param *sources, 
   phys_param *phys_param, 
   mtype **At)
{
  int k,m,s;
  int nnx, nny, nnz;
  dtype xmin, xmax, ymin, ymax, zmin, zmax, dx, dy, dz, A;
  dtype v;
  dtype greenr, greeni, phir, phii;
  dtype temp1r, temp1i, temp2r, temp2i;
  dtype phi_green_r, phi_green_i;
  dtype j_phi_green_r, j_phi_green_i;

  dtype tau_factor_r, tau_factor_i;
  dtype fre_gamma_r, fre_gamma_i;

  LINK2D sparse, src_element;
  LINK   det_list;

  nnx=phys_param->Ni;
  nny=phys_param->Nj;
  nnz=phys_param->Nl;

  v=phys_param->v;
  xmin=phys_param->xmin;
  xmax=phys_param->xmax;
  ymin=phys_param->ymin;
  ymax=phys_param->ymax;
  zmin=phys_param->zmin;
  zmax=phys_param->zmax;
 
  dx=(xmax-xmin)/((dtype)(nnx-1));
  dy=(ymax-ymin)/((dtype)(nny-1));
  dz=(zmax-zmin)/((dtype)(nnz-1));
  A = dx * dy * dz;


  sparse=phys_param->sparse;
  s=0;
  for(k=0; k<phys_param->K; k++){
    src_element=find_element_2d(sparse, k);
    if(src_element==NULL){
      fprintf(stderr, "calc_frechet_col: problem with linked list. Aborting.\n");
      return -1;
    }
    det_list=src_element->d;

    for(m=0; m<phys_param->M; m++) {
      if(find_data_in_list(det_list,m) >= 0){
        greenr=green[m][i][j][l][0];
        greeni=green[m][i][j][l][1];

        phir=phi_x[k][i][j][l][0];
        phii=phi_x[k][i][j][l][1];

       
        phi_green_r=cplxmultr(phir,phii, greenr, greeni);
        phi_green_i=cplxmulti(phir,phii, greenr, greeni);

        j_phi_green_r=-phi_green_i;
        j_phi_green_i=phi_green_r;

        tau_factor_r=1.0;
        tau_factor_i=-sources[0].omega*fluo[0][i][j][l]; 

        fre_gamma_r=cplxmultr(phi_green_r*A, phi_green_i*A, tau_factor_r, tau_factor_i );
        fre_gamma_i=cplxmulti(phi_green_r*A, phi_green_i*A, tau_factor_r, tau_factor_i );

        if(u==0) {
            At[s][0]=-j_phi_green_r*fluo[1][i][j][l]*sources[0].omega*A;
            At[s][1]=-j_phi_green_i*fluo[1][i][j][l]*sources[0].omega*A;
	}
        else if(u==1) {
            At[s][0]=fre_gamma_r;
            At[s][1]=fre_gamma_i; 
        }
	
        /*At[s][0]=-j_phi_green_r*fluo[1][i][j][l]*sources[0].omega*A;
        At[s][1]=-j_phi_green_i*fluo[1][i][j][l]*sources[0].omega*A; */

        s++;

      }
    }/* End for m */
  }/* End for k*/

  return 0;
}


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
)
{
  int k,m,s;
  int nnx, nny, nnz;
  dtype xmin, xmax, ymin, ymax, zmin, zmax, dx, dy, dz, A;
  dtype v;
  dtype greenr, greeni, phir, phii;
  dtype temp1r, temp1i, temp2r, temp2i;
  dtype phi_green_r, phi_green_i;
  dtype j_phi_green_r, j_phi_green_i;

  dtype tau_factor_r, tau_factor_i;
  dtype fre_gamma_r, fre_gamma_i;

  LINK2D sparse, src_element;
  LINK   det_list;

  nnx=phys_param->Ni;
  nny=phys_param->Nj;
  nnz=phys_param->Nl;

  v=phys_param->v;
  xmin=phys_param->xmin;
  xmax=phys_param->xmax;
  ymin=phys_param->ymin;
  ymax=phys_param->ymax;
  zmin=phys_param->zmin;
  zmax=phys_param->zmax;
 
  dx=(xmax-xmin)/((dtype)(nnx-1));
  dy=(ymax-ymin)/((dtype)(nny-1));
  dz=(zmax-zmin)/((dtype)(nnz-1));
  A = dx * dy * dz;


  sparse=phys_param->sparse;
  s=0;
  for(k=0; k<phys_param->K; k++){
    src_element=find_element_2d(sparse, k);
    if(src_element==NULL){
      fprintf(stderr, "calc_frechet_col: problem with linked list. Aborting.\n");
      return -1;
    }
    det_list=src_element->d;

    for(m=0; m<phys_param->M; m++) {
      if(find_data_in_list(det_list,m) >= 0){
        greenr=green[m][i][j][l][0];
        greeni=green[m][i][j][l][1];

        phir=phi_x[k][i][j][l][0];
        phii=phi_x[k][i][j][l][1];


        phi_green_r=cplxmultr(phir,phii, greenr, greeni); 
        phi_green_i=cplxmulti(phir,phii, greenr, greeni); 


        if(u==0) {
            At[s][0]=phi_green_r*A;
            At[s][1]=phi_green_i*A; 
        }
        else if(u==1) {
           j_phi_green_r=-phi_green_i;
           j_phi_green_i= phi_green_r;
            At[s][0]=j_phi_green_r*A;
            At[s][1]=j_phi_green_i*A;
        }

        s++;

      }
    }/* End for m */
  }/* End for k*/

  return 0;
}

int mult_frechet_row(
  /* Input */
  mtype *****phi,
  mtype *****green,
  mtype ****x,
  src_param *srcparam,
  det_param *detparam,
  phys_param *physparam
  /* Output */
){


   mtype ****At1;         /* At[u,i,j,l,c] */
   mtype ****At2;         /* At[u,i,j,l,c] */
   mtype ***Q;         /* At[u,i,j,l,c] */
   mtype At1r,At1i,At2i,At2r;
   FILE *fp=NULL;
  int nnx, nny, nnz,S,s1,s2;
  int i,j,l;
  int *k_s=NULL, *m_s=NULL;
 printf("Mult\n"); 
  nnx=physparam->Ni;
  nny=physparam->Nj;
  nnz=physparam->Nl;
  S=physparam->S;
  
  At1 = multialloc(sizeof(mtype), 4, nnx,nny,nnz,2);
  At2 = multialloc(sizeof(mtype), 4, nnx,nny,nnz,2);
  Q = multialloc(sizeof(mtype), 3, S, S,2);
  k_s=(int *)malloc(S*sizeof(int)); 
  m_s=(int *)malloc(S*sizeof(int)); 
 
  for(s1=0;s1<S; s1++){
    k_s[s1] = findk(physparam, s1);
    m_s[s1] = findm(physparam, s1);
  }

  for(s1=0;s1<S; s1++){
    printf("%d\n", s1);
    calc_frechet_row1(k_s[s1],m_s[s1],phi,green,x,srcparam,detparam,physparam,At1);
    for(s2=0;s2<S; s2++){
      Q[s1][s2][0]=0.;
      Q[s1][s2][1]=0.;
      calc_frechet_row1(k_s[s2],m_s[s2],phi,green,x,srcparam,detparam,physparam,At2);

      for (i=1; i<nnx-1; i++)
      for (j=1; j<nny-1; j++)
      for (l=1; l<nnz-1; l++) {
        At1r=At1[i][j][l][0];
        At1i=At1[i][j][l][1];
        At2r=At2[i][j][l][0];
        At2i=-At2[i][j][l][1];

        Q[s1][s2][0] += (At1r * At2r - At1i * At2i) ;
        Q[s1][s2][1] += (At1r * At2i + At1i * At2r) ;
      }
      
    }
  }

  fp = datOpen("Q.dat", "w+b");
  write_float_array(fp, "Q", &Q[0][0][0], 3, S, S, 2);
  datClose(fp);


  multifree(At1,4);
  multifree(At2,4);
  multifree(Q,3);
  free(k_s);
  free(m_s);

  return 0;
}

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
)
{
  int u,i,j,l;
  int nnx, nny, nnz;
  dtype xmin, xmax, ymin, ymax, zmin, zmax, dx, dy, dz, A;
  dtype greenr, greeni, phir, phii;
  double dgr, dgi, dpr, dpi;
 
  nnx=physparam->Ni;
  nny=physparam->Nj;
  nnz=physparam->Nl;
 
  xmin=physparam->xmin;
  xmax=physparam->xmax;
  ymin=physparam->ymin;
  ymax=physparam->ymax;
  zmin=physparam->zmin;
  zmax=physparam->zmax;
 
  dx=(xmax-xmin)/((dtype)(nnx-1));
  dy=(ymax-ymin)/((dtype)(nny-1));
  dz=(zmax-zmin)/((dtype)(nnz-1));
  A = dx * dy * dz;
 
 
/*  u = MUA;    */
  for (i=1; i<nnx-1; i++)
  for (j=1; j<nny-1; j++)
  for (l=1; l<nnz-1; l++) {
    greenr=green[m][i][j][l][0];
    greeni=green[m][i][j][l][1];
    phir=phi[k][i][j][l][0];
    phii=phi[k][i][j][l][1];
    At[i][j][l][0] =  (greenr * phir - greeni * phii) * 1; /* A */
    At[i][j][l][1] =  (greenr * phii + greeni * phir) * 1;
  }
 
/* 
  u = MUS;
  for (i=1; i<nnx-1; i++)
  for (j=1; j<nny-1; j++)
  for (l=1; l<nnz-1; l++) {
    At[u][i][j][l][0] = 0.0;
    At[u][i][j][l][1] = 0.0;
 
    dgr = (green[m][i][j][l][0] - green[m][i-1][j][l][0])/dx;
    dgi = (green[m][i][j][l][1] - green[m][i-1][j][l][1])/dx;
    dpr = (phi[k][i][j][l][0] - phi[k][i-1][j][l][0])/dx;
    dpi = (phi[k][i][j][l][1] - phi[k][i-1][j][l][1])/dx;
    At[u][i][j][l][0] += dgr * dpr - dgi * dpi;
    At[u][i][j][l][1] += dgi * dpr + dgr * dpi;
 
    dgr = (green[m][i+1][j][l][0] - green[m][i][j][l][0])/dx;
    dgi = (green[m][i+1][j][l][1] - green[m][i][j][l][1])/dx;
    dpr = (phi[k][i+1][j][l][0] - phi[k][i][j][l][0])/dx;
    dpi = (phi[k][i+1][j][l][1] - phi[k][i][j][l][1])/dx;
    At[u][i][j][l][0] += dgr * dpr - dgi * dpi;
    At[u][i][j][l][1] += dgi * dpr + dgr * dpi;
 
    dgr = (green[m][i][j][l][0] - green[m][i][j-1][l][0])/dy;
    dgi = (green[m][i][j][l][1] - green[m][i][j-1][l][1])/dy;
    dpr = (phi[k][i][j][l][0] - phi[k][i][j-1][l][0])/dy;
    dpi = (phi[k][i][j][l][1] - phi[k][i][j-1][l][1])/dy;
    At[u][i][j][l][0] += dgr * dpr - dgi * dpi;            

    dgr = (green[m][i][j+1][l][0] - green[m][i][j][l][0])/dy;
    dgi = (green[m][i][j+1][l][1] - green[m][i][j][l][1])/dy;
    dpr = (phi[k][i][j+1][l][0] - phi[k][i][j][l][0])/dy;
    dpi = (phi[k][i][j+1][l][1] - phi[k][i][j][l][1])/dy;
    At[u][i][j][l][0] += dgr * dpr - dgi * dpi;
    At[u][i][j][l][1] += dgi * dpr + dgr * dpi;
 
    dgr = (green[m][i][j][l][0] - green[m][i][j][l-1][0])/dz;
    dgi = (green[m][i][j][l][1] - green[m][i][j][l-1][1])/dz;
    dpr = (phi[k][i][j][l][0] - phi[k][i][j][l-1][0])/dz;
    dpi = (phi[k][i][j][l][1] - phi[k][i][j][l-1][1])/dz;
    At[u][i][j][l][0] += dgr * dpr - dgi * dpi;
    At[u][i][j][l][1] += dgi * dpr + dgr * dpi;
 
    dgr = (green[m][i][j][l+1][0] - green[m][i][j][l][0])/dz;
    dgi = (green[m][i][j][l+1][1] - green[m][i][j][l][1])/dz;
    dpr = (phi[k][i][j][l+1][0] - phi[k][i][j][l][0])/dz;
    dpi = (phi[k][i][j][l+1][1] - phi[k][i][j][l][1])/dz;
 
    At[u][i][j][l][0] += dgr * dpr - dgi * dpi;
    At[u][i][j][l][1] += dgi * dpr + dgr * dpi;
 
    At[u][i][j][l][0] *= -A / 2.0;
    At[u][i][j][l][1] *= -A / 2.0;
  }
*/ 
  return 0;
}                





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
)
{
  int k,m,u,i,j,l;
  int nnx, nny, nnz;
  dtype xmin, xmax, ymin, ymax, zmin, zmax, dx, dy, dz, A;
  dtype greenr, greeni, phir, phii;
  double dgr, dgi, dpr, dpi;
 
  nnx=physparam->Ni;
  nny=physparam->Nj;
  nnz=physparam->Nl;
 
  xmin=physparam->xmin;
  xmax=physparam->xmax;
  ymin=physparam->ymin;
  ymax=physparam->ymax;
  zmin=physparam->zmin;
  zmax=physparam->zmax;
 
  dx=(xmax-xmin)/((dtype)(nnx-1));
  dy=(ymax-ymin)/((dtype)(nny-1));
  dz=(zmax-zmin)/((dtype)(nnz-1));
  A = dx * dy * dz;
 
  k = findk(physparam, s);
  m = findm(physparam, s);
 
  u = MUA;    
  for (i=1; i<nnx-1; i++)
  for (j=1; j<nny-1; j++)
  for (l=1; l<nnz-1; l++) {
    greenr=green[m][i][j][l][0];
    greeni=green[m][i][j][l][1];
    phir=phi[k][i][j][l][0];
    phii=phi[k][i][j][l][1];
    At[u][i][j][l][0] = - (greenr * phir - greeni * phii) * A;
    At[u][i][j][l][1] = - (greenr * phii + greeni * phir) * A;
  }
 
  u = MUS;
 
  for (i=1; i<nnx-1; i++)
  for (j=1; j<nny-1; j++)
  for (l=1; l<nnz-1; l++) {
    At[u][i][j][l][0] = 0.0;
    At[u][i][j][l][1] = 0.0;
 
    dgr = (green[m][i][j][l][0] - green[m][i-1][j][l][0])/dx;
    dgi = (green[m][i][j][l][1] - green[m][i-1][j][l][1])/dx;
    dpr = (phi[k][i][j][l][0] - phi[k][i-1][j][l][0])/dx;
    dpi = (phi[k][i][j][l][1] - phi[k][i-1][j][l][1])/dx;
    At[u][i][j][l][0] += dgr * dpr - dgi * dpi;
    At[u][i][j][l][1] += dgi * dpr + dgr * dpi;
 
    dgr = (green[m][i+1][j][l][0] - green[m][i][j][l][0])/dx;
    dgi = (green[m][i+1][j][l][1] - green[m][i][j][l][1])/dx;
    dpr = (phi[k][i+1][j][l][0] - phi[k][i][j][l][0])/dx;
    dpi = (phi[k][i+1][j][l][1] - phi[k][i][j][l][1])/dx;
    At[u][i][j][l][0] += dgr * dpr - dgi * dpi;
    At[u][i][j][l][1] += dgi * dpr + dgr * dpi;
 
    dgr = (green[m][i][j][l][0] - green[m][i][j-1][l][0])/dy;
    dgi = (green[m][i][j][l][1] - green[m][i][j-1][l][1])/dy;
    dpr = (phi[k][i][j][l][0] - phi[k][i][j-1][l][0])/dy;
    dpi = (phi[k][i][j][l][1] - phi[k][i][j-1][l][1])/dy;
    At[u][i][j][l][0] += dgr * dpr - dgi * dpi;            

    dgr = (green[m][i][j+1][l][0] - green[m][i][j][l][0])/dy;
    dgi = (green[m][i][j+1][l][1] - green[m][i][j][l][1])/dy;
    dpr = (phi[k][i][j+1][l][0] - phi[k][i][j][l][0])/dy;
    dpi = (phi[k][i][j+1][l][1] - phi[k][i][j][l][1])/dy;
    At[u][i][j][l][0] += dgr * dpr - dgi * dpi;
    At[u][i][j][l][1] += dgi * dpr + dgr * dpi;
 
    dgr = (green[m][i][j][l][0] - green[m][i][j][l-1][0])/dz;
    dgi = (green[m][i][j][l][1] - green[m][i][j][l-1][1])/dz;
    dpr = (phi[k][i][j][l][0] - phi[k][i][j][l-1][0])/dz;
    dpi = (phi[k][i][j][l][1] - phi[k][i][j][l-1][1])/dz;
    At[u][i][j][l][0] += dgr * dpr - dgi * dpi;
    At[u][i][j][l][1] += dgi * dpr + dgr * dpi;
 
    dgr = (green[m][i][j][l+1][0] - green[m][i][j][l][0])/dz;
    dgi = (green[m][i][j][l+1][1] - green[m][i][j][l][1])/dz;
    dpr = (phi[k][i][j][l+1][0] - phi[k][i][j][l][0])/dz;
    dpi = (phi[k][i][j][l+1][1] - phi[k][i][j][l][1])/dz;
 
    At[u][i][j][l][0] += dgr * dpr - dgi * dpi;
    At[u][i][j][l][1] += dgi * dpr + dgr * dpi;
 
    At[u][i][j][l][0] *= -A / 2.0;
    At[u][i][j][l][1] *= -A / 2.0;
  }
 
  return 0;
}                

