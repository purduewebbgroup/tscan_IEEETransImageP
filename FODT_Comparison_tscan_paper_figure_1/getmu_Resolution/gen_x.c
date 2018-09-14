#include "defs.h"

int main(int argc, char *argv[]){
  char datname[255] ; 
  char varname[255] ; 
  FILE *fpd;
  mtype mu[2][65][65][65];  // modify: Ni where i = x,y,z
  //double x, y, z,r;
  int i, j, l;
  

  if(argc!=2){
    printf("Usage: %s filename\n",argv[0]);
    exit(-1);
  }

  strcpy(datname, argv[1]);
  strcpy(varname, "mu");
  strcat(datname,".dat");

  for(i=0; i<65; i++) // modify: Ni where i = x,y,z
  for(j=0; j<65; j++)
  for(l=0; l<65; l++)
  {	 
     mu[0][i][j][l]=0.0165;   // D;   alpha = 1
     mu[1][i][j][l]=0.2;      // mua; alpha = 1
  }

  fpd = datOpen(datname, "w+b");
  write_float_array(fpd, varname, &mu[0][0][0][0], 4, 2, 65,65,65); // modify: Ni where i = x,y,z
  datClose(fpd);

  return 0;
}


