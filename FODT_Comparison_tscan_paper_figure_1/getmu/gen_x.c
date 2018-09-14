#include "defs.h"

int main(int argc, char *argv[]){
  char datname[255] ; 
  char varname[255] ; 
  FILE *fpd;
  mtype mu[2][33][33][17]; 
  double x, y, z,r;
  int i, j, l;
  

  if(argc!=2){
    printf("Usage: %s filename\n",argv[0]);
    exit(-1);
  }

  strcpy(datname, argv[1]);
  strcpy(varname, "mu");
  strcat(datname,".dat");

  for(i=0; i<33; i++)
  for(j=0; j<33; j++)
  for(l=0; l<17; l++)
  {

     x=-8.18+(16.36/32.)*((double)i);
     y=-8.18+(16.36/32.)*((double)j);
     z=-3.03+(6.06/16.)*((double)l);
     mu[0][i][j][l]=0.047; 
     mu[1][i][j][l]=0.027;
     r=sqrt(x*x+y*y+z*z);
     if(r<=0.5){
     }
  }

  fpd = datOpen(datname, "w+b");
  write_float_array(fpd, varname, &mu[0][0][0][0], 4, 2,33,33,17);
  datClose(fpd);

  return 0;
}


