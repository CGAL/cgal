/*********************************************************/ 
/* TAUCS						 */ 
/* Author: Sivan Toledo					 */ 
/*********************************************************/ 

#include <stdio.h> 

int main(int argc, char* argv[]) 
{ 
  double A = 2.0; 
  double B = 3.0; 
  double C = 5.0; 
  double alpha = 7.0;
  double beta  = 11.0;
  int n = 1; 
  int ld = 1; 
  FILE * file;

  void dgemm (char*,char*,int*,
	      int*,int*,
	      double*,
	      double*,int*,double*,int*,
	      double*,
	      double*,int*);
  
  dgemm ("N","N", &n,&n,&n, &alpha, &A,&ld, &B,&ld, &beta, &C,&ld);

  printf("\n\n");
  printf("Linking with dgemm succedded\n");
  printf("\n\n");

  if (C == 97.0) {
    file = fopen( argv[1], "a" );
    if ( file == NULL ) {
      printf("Problem opening file.\n");
      return -1;
    }
    fprintf(file, "/* Definition for BLAS functions */\n");
    fprintf(file, "#define TAUCS_BLAS_NOUNDERSCORE\n");
    fclose( file );
    return 0;
  } else {
    return 1;
  }
} 
