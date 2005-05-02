/*********************************************************/ 
/* TAUCS						 */ 
/* Author: Sivan Toledo					 */ 
/*********************************************************/ 

#include <stdio.h> 
#include <stdlib.h> 
#include <math.h> 
#include <complex.h> 

int main(int argc, char* argv[]) 
{ 
  _Complex double x = -1.0;
  _Complex double y;
  FILE * file;
  
  y = csqrt(x) + x;

  /* we now test whether the real and imaginary parts can serve as lvalues */

  printf("\n\n");
  printf("C99 complex numbers seem to be supported, 1+sqrt(-1)=%f+%fi\n",
	 creal(y),cimag(y));
  printf("\n\n");

  file = fopen( argv[1], "a" );
  if ( file == NULL ) {
    printf("Problem opening file.\n");
    return 1;
  }
  fprintf(file, "/* Does the compiler support C99 complex numbers? */\n");
  fprintf(file, "#define TAUCS_C99_COMPLEX\n");
  fclose( file );

  return 0;
} 
