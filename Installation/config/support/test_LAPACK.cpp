// Test if LAPACK is available

#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#ifdef CGAL_USE_F2C
#define my_dgelss dgelss_
#else  
#define my_dgelss dgelss
#endif

extern "C" {
int my_dgelss(int *m, int *n, int *nrhs, 
		    double *a, int *lda, double *b, int *ldb, double *
		    s, double *rcond, int *rank, double *work, int *lwork,
		    int *info); 
 }


int main(int argc, char* argv[])
{
  double M=1,
    B=1;
  int m = 1,
    n = 1,
    nrhs = 1,
    lda = m,
    ldb = m,
    rank,
    lwork = 5*m,
    info;
  //c style
  double* sing_values = (double*) malloc(n*sizeof(double));
  double* work = (double*) malloc(lwork*sizeof(double));

  double rcond = -1;

  my_dgelss(&m, &n, &nrhs, &M, &lda, &B, &ldb, sing_values, 
	  &rcond, &rank, work, &lwork, &info);
  assert(info==0);
  assert(B==1.);
  //clean up 
  free(sing_values);
  free(work);

  std::cout << "ok for lapack" << std::endl;
  
  return 0;
}
