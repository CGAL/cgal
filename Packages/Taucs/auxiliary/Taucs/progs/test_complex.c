/*********************************************************/ 
/* TAUCS						 */ 
/* Author: Sivan Toledo					 */ 
/*********************************************************/ 

/* 

TAUCS_CONFIG_BEGIN
TAUCS_CONFIG_DEFAULT OFF
TAUCS_CONFIG TIMING
TAUCS_CONFIG BASE
TAUCS_CONFIG DREAL
TAUCS_CONFIG SREAL
TAUCS_CONFIG DCOMPLEX
TAUCS_CONFIG SCOMPLEX
TAUCS_CONFIG MALLOC_STUBS
TAUCS_CONFIG AD_HOC_TEST
TAUCS_CONFIG_END

*/

#include <stdio.h> 
#include <stdlib.h> 
#include <math.h> 
#include <complex.h> 

#include <taucs.h> 

int main(int argc, char* argv[]) 
{ 
  taucs_dcomplex m1 = taucs_zminusone;
  taucs_dcomplex p1 = taucs_zone;
  taucs_dcomplex z  = taucs_zzero;
  taucs_dcomplex i  = taucs_zsqrt(m1);
  taucs_dcomplex x  = taucs_zadd(p1,i);
  
  printf("\n\n");
#ifdef TAUCS_C99_COMPLEX
  printf("This build of TAUCS uses C99 complex-number support\n");
#else
  printf("This build of TAUCS uses TAUCS's built-in complex-number support\n");
#endif

  printf("\n");

  printf("+1       = %f+%fi\n",taucs_zreal(p1),taucs_zimag(p1));
  printf("-1       = %f+%fi\n",taucs_zreal(m1),taucs_zimag(m1));
  printf("0        = %f+%fi\n",taucs_zreal(z),taucs_zimag(z));
  printf("i        = %f+%fi\n",taucs_zreal(i),taucs_zimag(i));
  printf("1+i      = %f+%fi\n",taucs_zreal(x),taucs_zimag(x));

  {
    unsigned int   j;
    double         q = 1.0;
    taucs_dcomplex p = x;
    double         t1,t2;
    
    t1=taucs_wtime();
    for (j=0; j<10000000; j++) q+=1.0;
    t2=taucs_wtime();
    printf("real    addition performance approx %.2e Mflop/s\n",1e7/(t2-t1));

    t1=taucs_wtime();
    for (j=0; j<10000000; j++) p=taucs_zadd(p,x);
    t2=taucs_wtime();
    printf("complex addition performance approx %.2e Mflop/s\n",2e7/(t2-t1));
  }

  printf("\n\n");

  return 0;
} 
