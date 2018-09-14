#include "defs.h"

int main(int argc, char *argv[]){
  char datname[255] ; 
  char varname[255] ; 
  FILE *fpd;
  mtype mu[2][33][33][17]; 
  double x, y, z,r1,r2;
  int i, j, l;
  

  if(argc!=2){
    printf("Usage: %s filename\n",argv[0]);
    exit(-1);
  }

  strcpy(datname, argv[1]);
  strcpy(varname, "fluo");
  strcat(datname,".dat");

  for(i=0; i<33; i++)
  for(j=0; j<33; j++)
  for(l=0; l<17; l++)
  {

     x=-4.18+(8.36/32.)*((double)i);
     y=-4.18+(8.36/32.)*((double)j);
     z=-3.03+(6.06/16.)*((double)l);
     mu[0][i][j][l]=0.0; 
     mu[1][i][j][l]=0.0;
     r1=sqrt((x-1)*(x-1)+(y-1.0)*(y-1.0)+(z-1.1)*(z-1.1));
     r2=sqrt((x+1.30)*(x+1.30)+(y+1.1)*(y+1.1)+(z+.9)*(z+.9));
     if(r1<=.7012){
        mu[0][i][j][l]=0.55e-9; 
        mu[1][i][j][l]= 1;
     }
     if(r2<=.7012){
        mu[0][i][j][l]=0.55e-9; 
        mu[1][i][j][l]= 1;
     }
  }

  fpd = datOpen(datname, "w+b");
  write_float_array(fpd, varname, &mu[0][0][0][0], 4, 2,33,33,17);
  datClose(fpd);

  return 0;
}


