/*********************************************************/ 
/* TAUCS						 */ 
/* Author: Sivan Toledo					 */ 
/*********************************************************/ 

#include <stdio.h> 
#include <stdlib.h> 
#include <cilk.h> 

#pragma lang -C

cilk int f()
{
  return 0;
}

cilk int main(int argc, char* argv[]) 
{
  int x;
  FILE* file;

  x = spawn f();

  sync;

  printf("\n\n");
  printf("Compiling Cilk succedded\n");
  printf("\n\n");

  file = fopen( argv[1], "a" );
  if ( file == NULL ) {
    printf("Problem opening file.\n");
    return 1;
  }

  fprintf(file, "/* Do we have a Cilk compiler? */\n");
  fprintf(file, "#define TAUCS_CILK\n");
  fclose( file );

  return 0;
}
