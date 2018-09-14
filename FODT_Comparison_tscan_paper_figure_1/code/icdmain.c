/*
 * ICD-Born iterative method for optical diffusion tomography 
 * Jong's paper, JOSA-A, Oct, 1999
 * 
 * Seungseok Oh 
 * 10/25/2000
 */

#include <stdio.h>
#include <math.h>
#include "allocate.h"
#include "defs.h"


int main(int argc, char **argv){
  int i,j,l,c,s,k,m;
  int status;
  int fwd_flag=0;
  phys_param  *physparam;
  src_param    *srcparam;
  det_param    *detparam;

  prior_param *priorparam_tau;
  prior_param *priorparam_gamma;

  prior_param *priorparam_D;
  prior_param *priorparam_mua;

  dtype *x, *y, *z, dx, dy, dz;
  dtype alpha_fixed=9.4e-6;
  dtype max_snr=100.0; 
  mtype ****mu_x, ****mu_m, ****fluo, ****fluohat;  
  mtype **meas, **meas_true, **yerror;
  mtype **meas_m;
  mtype *****phi_x;
  mtype *****phi_m;
  mtype *snr;
  
  FILE *fp=NULL;
  char physfile[255];

  char priorfile_D[255];
  char priorfile_mua[255]; 

  char priorfile_tau[255]; 
  char priorfile_gamma[255];
  char measfile[255];
  char phantomfile_x[255];
  char phantomfile_m[255];
  char fluofile[255];
/*  config_param config; */
  fluo_config_param config;
  double tempr, tempi, temp2r, temp2i, temp3r, temp3i;    
  dtype b[26]; 

  //max_snr=1000000000000000000000.0;
  max_snr=30.0;

  if (argc!=2 && argc!=9) {
    fprintf(stderr, "Usage: %s configfile \n", argv[0]);
    fprintf(stderr, "OR   : %s -F physfile phantom_ex phantom_em fluor alpha wgtr wgti\n", argv[0]);
/*  max_snr sets an upper limit on the SNR for the generated noisy data.  Here, I am disabling this from the 
     argument list and setting a default value to reduce confusion for the new user */
    exit(1);
  }

  if(argc==2){ 
    fp = fopen(argv[1],"r");
    if (fp==NULL) {
      fprintf(stderr, "ERROR! Cannot load configuration file!\n");
      exit(1);
    }
    get_config_fluo(fp, &config);
    fclose(fp);

/*    strcpy(physfile, config.phys_file); */
    strcpy(priorfile_D, config.prior_D_file);
    strcpy(priorfile_mua, config.prior_mua_file);
/*   strcpy(measfile, config.meas_file); */

    strcpy(physfile, config.phys_file);
    strcpy(priorfile_tau, config.prior_tau_file);
    strcpy(priorfile_gamma, config.prior_gamma_file);
    strcpy(measfile, config.meas_file); 
    strcpy(phantomfile_x, config.mu_x_file);
/*    strcpy(phantomfile_m, config.mu_m_file); */

  }
  else{
    fwd_flag=1;
    strcpy(physfile, argv[2]);
    strcpy(phantomfile_x, argv[3]);
    strcpy(phantomfile_m, argv[4]);
    strcpy(fluofile, argv[5]);
    alpha_fixed=atof(argv[6]);
  }

  printf("Reading in physical parameters from file `%s`...\n",physfile);
  fp=fopen(physfile, "r");
  if(fp==NULL){
    fprintf(stderr,"Error: can't open `%s`\n", physfile);
    exit(1);
  }
  
  physparam=alloc_phys_param();
  get_all_phys(fp, physparam, &srcparam, &detparam);  
  fclose(fp);

/* Milstein - BEGIN */

  if(!fwd_flag){
    physparam->wgtr=config.init_wgtr;
    physparam->wgti=0.0;
  }
  else{
    physparam->wgtr=atof(argv[7]);
    physparam->wgti=atof(argv[8]);
  }    

/* Milstein - END */   


  /* make a grid */
  dx=(physparam->xmax - physparam->xmin) / ((mtype)(physparam->Ni - 1));
  dy=(physparam->ymax - physparam->ymin) / ((mtype)(physparam->Nj - 1));
  dz=(physparam->zmax - physparam->zmin) / ((mtype)(physparam->Nl - 1));

  x=(dtype *)malloc((physparam->Ni)*sizeof(dtype));
  y=(dtype *)malloc((physparam->Nj)*sizeof(dtype));
  z=(dtype *)malloc((physparam->Nl)*sizeof(dtype));

  for(i=0; i<physparam->Ni; i++)
    x[i]= physparam->xmin + ((mtype)(i)) * dx;
  for(j=0; j<physparam->Nj; j++)
    y[j]= physparam->ymin + ((mtype)(j)) * dy;
  for(l=0; l<physparam->Nl; l++)
    z[l]= physparam->zmin + ((mtype)(l)) * dz;


  /* Set prior model parameters */

  if(!fwd_flag){
    priorparam_tau = alloc_prior_param(NULL, 0, 0, 26, b);    
    printf("Reading in model parameters from file `%s`...\n",priorfile_tau);
    fp=fopen(priorfile_tau,"r");
    if(fp==NULL){
      fprintf(stderr,"Error: can't open `%s`\n", priorfile_tau);
      exit(1);
    }
    get_prior(fp,priorparam_tau);
    fclose(fp);

    priorparam_gamma = alloc_prior_param(NULL, 0, 0, 26, b);    
    printf("Reading in model parameters from file `%s`...\n",priorfile_gamma);
    fp=fopen(priorfile_gamma,"r");
    if(fp==NULL){
      fprintf(stderr,"Error: can't open `%s`\n", priorfile_gamma);
      exit(1);
    }
    get_prior(fp,priorparam_gamma);
    fclose(fp);

    priorparam_D = alloc_prior_param(NULL, 0, 0, 26, b);    
    printf("Reading in model parameters from file `%s`...\n",priorfile_D);
    fp=fopen(priorfile_D,"r");
    if(fp==NULL){
      fprintf(stderr,"Error: can't open `%s`\n", priorfile_D);
      exit(1);
    }
    get_prior(fp,priorparam_D);
    fclose(fp);

/* This code has the ability to perform update iterations on mu_m and D_m, but
    it is not recommended!  It serves little purpose, and it has rarely been
    tested in a while */
    priorparam_mua = alloc_prior_param(NULL, 0, 0, 26, b);   
    printf("Reading in model parameters from file `%s`...\n",priorfile_mua);
    fp=fopen(priorfile_mua,"r");
    if(fp==NULL){
      fprintf(stderr,"Error: can't open `%s`\n", priorfile_mua);
      exit(1);
    }
    get_prior(fp,priorparam_mua);
    fclose(fp);


  }

  /* Set absorption and scattering coeff. */

  mu_x = multialloc(sizeof(mtype), 4, 2, physparam->Ni, physparam->Nj, physparam->Nl);
  mu_m = multialloc(sizeof(mtype), 4, 2, physparam->Ni, physparam->Nj, physparam->Nl);
  fluo = multialloc(sizeof(mtype), 4, 2, physparam->Ni, physparam->Nj, physparam->Nl);
  fluohat= multialloc(sizeof(mtype), 4, 2, physparam->Ni, physparam->Nj, physparam->Nl);

/****  if(fwd_flag){  */
  fp = datOpen(phantomfile_x, "r");
  if(fp==NULL){
    fprintf(stderr,"Error: can't open `%s`\n", phantomfile_x);
    exit(1);
  }
  status=read_float_array(fp, "mu" , &mu_x[0][0][0][0] , 4, 2, physparam->Ni, physparam->Nj, physparam->Nl );

  if(status==FOERR){
    fprintf(stderr,"Error reading mu_x. Aborting.\n");
    exit(1);
  }
  datClose(fp);

  if(fwd_flag){  
    fp = datOpen(phantomfile_m, "r");
    if(fp==NULL){
      fprintf(stderr,"Error: can't open `%s`\n", phantomfile_m);
      exit(1);
    }
    status=read_float_array(fp, "mu" , &mu_m[0][0][0][0] , 4, 2, physparam->Ni, physparam->Nj, physparam->Nl );

    if(status==FOERR){
      fprintf(stderr,"Error reading mu_m. Aborting.\n");
      exit(1);
    }
    datClose(fp);
    
    fp = datOpen(fluofile, "r");
    if(fp==NULL){
      fprintf(stderr,"Error: can't open `%s`\n", fluofile);
      exit(1);
    }
    status=read_float_array(fp, "fluo" , &fluo[0][0][0][0] , 4, 2, physparam->Ni, physparam->Nj, physparam->Nl );

    if(status==FOERR){
      fprintf(stderr,"Error reading fluorescent parameters. Aborting.\n");
      exit(1);
    }
    datClose(fp);

  }

  /* simulate measurement */

  phi_x    = multialloc(sizeof(mtype), 5, physparam->K, physparam->Ni, 
                      physparam->Nj, physparam->Nl, 2);

  phi_m  = multialloc(sizeof(mtype), 5, physparam->K, physparam->Ni, 
                      physparam->Nj, physparam->Nl, 2);
  
  meas   = multialloc(sizeof(mtype), 2, physparam->S, 2);
  meas_m = multialloc(sizeof(mtype), 2, physparam->S, 2);

  yerror = multialloc(sizeof(mtype), 2, physparam->S, 2);
  meas_true   = multialloc(sizeof(mtype), 2, physparam->S, 2);
  snr   = malloc(sizeof(mtype)*physparam->S);

  if(srcparam==NULL || detparam==NULL|| physparam==NULL|| 
     phi_x==NULL || meas==NULL || yerror==NULL || meas_true==NULL||
     meas_m==NULL || phi_m==NULL)  {
    printf("Memory problems.\n");
    exit(1);
  }

  if(fwd_flag){
     printf("Calling calc_phi...\n");
     calc_phi(srcparam, detparam, mu_x, physparam, phi_x, meas);
     printf("Calling calc_phi_fluo...\n");
     calc_phi_fluo(fluo, phi_x, srcparam, detparam, mu_m, physparam, phi_m, meas_m);

     printf("Adding detector noise...\n");
     add_detector_noise(physparam, alpha_fixed, max_snr, meas_m, snr);

     for (s=0; s<physparam->S; s++) {
       k = findk(physparam,s);
       m = findm(physparam,s);

       tempr = cplxmultr(srcparam[k].calir, srcparam[k].calii,
                         detparam[m].calir, detparam[m].calii);
       tempi = cplxmulti(srcparam[k].calir, srcparam[k].calii,
                         detparam[m].calir, detparam[m].calii);
       temp2r = cplxmultr(meas_m[s][0], meas_m[s][1], tempr, tempi);
       temp2i = cplxmulti(meas_m[s][0], meas_m[s][1], tempr, tempi);

       /* Milstein - Begin */
       temp3r = cplxmultr(temp2r, temp2i, physparam->wgtr, physparam->wgti);
       temp3i = cplxmulti(temp2r, temp2i, physparam->wgtr, physparam->wgti);
 
       fprintf(stderr, "%2d  %2d  %2d : %f  %f  %f  %f  %f  %f \n",
                 s, k, m, tempr, tempi, meas_m[s][0], meas_m[s][1], temp3r, temp3i);
       meas_m[s][0] = temp3r;
       meas_m[s][1] = temp3i;
 
       /* Milstein - END */    
     }

     fp = datOpen("meas.dat", "w+b");
     write_float_array(fp, "meas", &meas_m[0][0], 2, physparam->S, 2);
     datClose(fp);
     fp = datOpen("snr.dat", "w+b");
     write_float_array(fp, "snr", &snr[0], 1, physparam->S);
     datClose(fp);
  }

  else{
    printf("Reading measurement data...\n");
    fp = datOpen(measfile, "r+b");
    status=read_float_array(fp, "meas" , &meas_true[0][0] , 2, physparam->S, 2);
    if(status==FOERR){
      fprintf(stderr,"Error reading measurements. Aborting.\n");
      exit(1);
    }
    datClose(fp);

    /* initial guess of tau and gamma , and mua_m and D_m */ 

    for(i=0; i<physparam->Ni; i++)
    for(j=0; j<physparam->Nj; j++)
    for(l=0; l<physparam->Nl; l++) {   
      fluohat[0][i][j][l] = config.tau_backg; 
      fluohat[1][i][j][l] = config.etamu_backg; 
      mu_m[0][i][j][l] = config.mua_backg; 
      mu_m[1][i][j][l] = config.D_backg; 
    } 

    if(config.init_guess_fluo_flag){
      fp = datOpen(config.init_guess_fluo_path, "r");
      if(fp==NULL){
        fprintf(stderr,"Error: can't open `%s`\n", config.init_guess_fluo_path);
        exit(1);
      }
      status=read_float_array(fp, config.init_guess_fluo_varname , 
                              &fluohat[0][0][0][0] , 4, 2, 
                              physparam->Ni, physparam->Nj, physparam->Nl );
 
      if(status==FOERR){
        fprintf(stderr,"Error reading initial guess file. Aborting.\n");
        exit(1);
      }
      datClose(fp);
    }

    if(config.init_guess_opt_flag){
      fp = datOpen(config.init_guess_opt_path, "r");
      if(fp==NULL){
        fprintf(stderr,"Error: can't open `%s`\n", config.init_guess_opt_path);
        exit(1);
      }
      status=read_float_array(fp, config.init_guess_opt_varname , 
                              &mu_m[0][0][0][0] , 4, 2, 
                              physparam->Ni, physparam->Nj, physparam->Nl );
 
      if(status==FOERR){
        fprintf(stderr,"Error reading initial guess file. Aborting.\n");
        exit(1);
      }
      datClose(fp);
    }

    /* Intitialize the source-detector weights. 
       Initial guess of calibration is 1.0+0.0j. */
/*
    for (k=0; k<physparam->K; k++) {
      srcparam[k].calir = 1.0;
      srcparam[k].calii = 0.0;
    }

    for (m=0; m<physparam->M; m++) {
      detparam[m].calir = 1.0;
      detparam[m].calii = 0.0;
    }  
*/

    for (k=0; k<physparam->K; k++) {
      printf("%f,%f\n",srcparam[k].calir,srcparam[k].calii);
    }

    for (m=0; m<physparam->M; m++) {
      printf("%f,%f\n", detparam[m].calir,detparam[m].calii);
    }  

    printf("Forward solution...\n");

    calc_phi(srcparam, detparam, mu_x, physparam, phi_x, meas);
    calc_phi_fluo(fluohat, phi_x, srcparam, detparam, mu_m, physparam, phi_m, meas_m);

    for (s=0; s<physparam->S; s++) 
    for (c=0; c<2; c++) {
      yerror[s][c] = meas_true[s][c] - meas_m[s][c];
      printf("yerror[%d][%d]=%e\n",s,c,yerror[s][c]);
    }


/*
    for (s=0; s<physparam->S; s++) {
       k = findk(physparam,s);
       m = findm(physparam,s);

       tempr = cplxmultr(srcparam[k].calir, srcparam[k].calii,
                         detparam[m].calir, detparam[m].calii);
       tempi = cplxmulti(srcparam[k].calir, srcparam[k].calii,
                         detparam[m].calir, detparam[m].calii);
       temp2r = cplxmultr(meas_m[s][0], meas_m[s][1], tempr, tempi);
       temp2i = cplxmulti(meas_m[s][0], meas_m[s][1], tempr, tempi);

       temp3r = cplxmultr(temp2r, temp2i, physparam->wgtr, physparam->wgti);
       temp3i = cplxmulti(temp2r, temp2i, physparam->wgtr, physparam->wgti);
 

       yerror[s][0] = meas_true[s][0] - temp3r;
       yerror[s][1] = meas_true[s][1] - temp3i;
       printf("yerror[%d][%d]=%e\n",s,0,yerror[s][0]);
    }
*/
    /* perform ICD iterations */

    printf("Inversion...\n");
    rvssolver_sg_fluo(fluohat, mu_m, phi_x, phi_m, meas_m, 
                 meas_true, physparam, srcparam, detparam, 
                 priorparam_tau, priorparam_gamma, 
                 priorparam_D, priorparam_mua, 
                 yerror, &config); 

  }


  delete_phys_param(physparam); 
  free(srcparam);
  free(detparam);

  if(!fwd_flag) {
    free_prior_param(priorparam_tau);
    free_prior_param(priorparam_gamma);
  }

  multifree(mu_x,4);
  multifree(mu_m,4);
  multifree(fluo,4);
  multifree(fluohat,4);
  multifree(phi_x,5);
  multifree(phi_m,5);
  multifree(meas,2);
  multifree(meas_m,2);
  multifree(meas_true,2);
  multifree(yerror,2);
  free(x);
  free(y);
  free(z);
  
  return 0;
}
