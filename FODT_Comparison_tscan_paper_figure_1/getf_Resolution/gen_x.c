#include "defs.h"

int main(int argc, char *argv[]){
  char datname[255] ; 
  char varname[255] ; 
  FILE *fpd;
  mtype mu[2][65][65][65];  // modify: Ni where i = x,y,z
  double x, y, z, alpha, mus, l_mfp, tolerance, delta;
  double xmin, ymin, zmin, xmax, ymax, zmax, Bounded;
  double delta_x, delta_y, delta_z;
  int Ni, Nj, Nl, i, j, l;
  

  if(argc!=2){
    printf("Usage: %s filename\n",argv[0]);
    exit(-1);
  }

  strcpy(datname, argv[1]);
  strcpy(varname, "fluo");
  strcat(datname,".dat");
  
  Ni = 65;
  Nj = 65;
  Nl = 65;
  
  alpha = 1.0;  // scaling parameter
  mus = 20.0/alpha;
  l_mfp = 1.0/mus;
  
  Bounded = 1;
  xmin = -l_mfp*Bounded;
  xmax = 32*l_mfp + l_mfp*Bounded;
  ymin = -l_mfp*Bounded;
  ymax = 32*l_mfp + l_mfp*Bounded;
  zmin = -l_mfp*Bounded;
  zmax = 32*l_mfp + l_mfp*Bounded;
  
  delta_x = (xmax-xmin)/(Ni-1);
  delta_y = (ymax-ymin)/(Nj-1);
  delta_z = (zmax-zmin)/(Nl-1);
  
  printf("\ndelta_x = %f\n", delta_x);
  printf("delta_y = %f\n", delta_y);
  printf("delta_z = %f\n", delta_z);
  
  tolerance = delta_x/2;
  
  delta = 4.0;

  for(i=0; i<Ni; i++)
  for(j=0; j<Nj; j++)
  for(l=0; l<Nl; l++)
  {	 
	 x=xmin+delta_x*((double)i);
     y=ymin+delta_y*((double)j); 
     z=zmin+delta_z*((double)l);
	 
     mu[0][i][j][l]=0.0; 
     mu[1][i][j][l]=0.0;
	 
	 // turn on four fluorescent voxels (point sources)
	 if(fabs(x-15.0*l_mfp)<tolerance && fabs(y-16.0*l_mfp)<tolerance && fabs(z-15.0*l_mfp)<tolerance){
		mu[0][i][j][l]=0.5e-9;  // Tau
        mu[1][i][j][l]= 1.0;      // Eta
		
		printf("\nx = %f\n", x);
		printf("y = %f\n", y);
		printf("z = %f\n", z);
	 }
	 
	 if(fabs(x-15.0*l_mfp)<tolerance && fabs(y-16.0*l_mfp)<tolerance && fabs(z-15.0*l_mfp-delta*delta_z)<tolerance){
		mu[0][i][j][l]=0.5e-9;  // Tau
        mu[1][i][j][l]= 1.5;      // Eta
		
		printf("\nx = %f\n", x);
		printf("y = %f\n", y);
		printf("z = %f\n", z);
	 }
	 
	 if(fabs(x-15.0*l_mfp-delta*delta_x)<tolerance && fabs(y-16.0*l_mfp)<tolerance && fabs(z-15.0*l_mfp)<tolerance){
		mu[0][i][j][l]=0.5e-9;  // Tau
        mu[1][i][j][l]= 0.75;      // Eta
		
		printf("\nx = %f\n", x);
		printf("y = %f\n", y);
		printf("z = %f\n", z);
	 }
	 
	 if(fabs(x-15.0*l_mfp-delta*delta_x)<tolerance && fabs(y-16.0*l_mfp)<tolerance && fabs(z-15.0*l_mfp-delta*delta_z)<tolerance){
		mu[0][i][j][l]=0.5e-9;  // Tau
        mu[1][i][j][l]= 0.5;      // Eta
		
		printf("\nx = %f\n", x);
		printf("y = %f\n", y);
		printf("z = %f\n", z);
	 }
	 

	 
  }

  fpd = datOpen(datname, "w+b");
  write_float_array(fpd, varname, &mu[0][0][0][0], 4, 2, Ni,Nj,Nl);
  datClose(fpd);

  return 0;
}
