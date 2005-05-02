/*********************************************************/ 
/*				TAUCS									 */ 
/*			Author: Sivan Toledo						 */ 
/*********************************************************/ 

#include <stdio.h> 

#define TAUCS_CORE_DOUBLE
#define TAUCS_IGNORE_CONFIG_TESTS

#include <taucs.h>

/*********************************************************/ 
/*														 */ 
/*********************************************************/ 


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
  
#ifdef TAUCS_BLAS_UNDERSCORE
  char * str = "#define TAUCS_BLAS_UNDERSCORE\n";
#else
  char * str = "#ifdef TAUCS_BLAS_UNDERSCORE\n#undef TAUCS_BLAS_UNDERSCORE\n#endif\n";
#endif
	
  taucs_gemm("N","N", &n,&n,&n, &alpha, &A,&ld, &B,&ld, &beta, &C,&ld);

  if (C == 97.0) {
    file = fopen( argv[1], "w" );
    if ( file == NULL ) {
      printf("Problem opening file.\n");
      return -1;
    }
    fprintf(file, "/* Definition for BLAS functions */\n");
    fprintf(file, "#define TAUCS_TEST_BLAS_UNDERSCORE\n");
    fprintf(file, "%s", str );
    fclose( file );
    return 0;
  } else {
    return 1;
  }
} 
