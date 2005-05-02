/*********************************************************/
/* TAUCS                                                 */
/* Author: Sivan Toledo                                  */
/*********************************************************/
/* 

TAUCS_CONFIG_BEGIN
TAUCS_CONFIG_DEFAULT OFF
TAUCS_CONFIG BASE
xxxTAUCS_CONFIG CILK
TAUCS_CONFIG DREAL
TAUCS_CONFIG FACTOR
TAUCS_CONFIG OOC_LLT
TAUCS_CONFIG METIS
TAUCS_CONFIG GENMMD
TAUCS_CONFIG COLAMD
TAUCS_CONFIG AMD
TAUCS_CONFIG GENERIC_COMPLEX
TAUCS_CONFIG MALLOC_STUBS
TAUCS_CONFIG MATRIX_GENERATORS
TAUCS_CONFIG AD_HOC_TEST
TAUCS_CONFIG_END

*/

#include <stdio.h>
#include <stdlib.h>
#include "taucs.h"

#ifdef TAUCS_BLAS_UNDERSCORE
#define my_dnrm2 dnrm2_
#else
#define my_dnrm2 dnrm2
#endif

double my_dnrm2();

int rnorm(taucs_ccs_matrix* A, double* x, double* y, double* aux)
{
  int i;
  int one = 1;
  double relerr;
  
  for(i=0; i<A->n; i++) aux[i] = y[i]-x[i];
  relerr = my_dnrm2(&(A->n),aux,&one)/my_dnrm2(&(A->n),x,&one); 
  taucs_printf("2-norm relative forward error %.2e \n",relerr);
  
  if (relerr > 1e-8) {
    taucs_printf("test failed, error seems too high (but\n");
    taucs_printf("maybe the matrix is too ill-conditioned)\n");
    return 1;
  }
  return 0;
}

int test_spd_orderings(taucs_ccs_matrix* A,
		       double* x, double* y, double* b, double* z)
{
  int rc;
  char* metis[] = {"taucs.factor.ordering=metis", NULL};
  char* genmmd[] = {"taucs.factor.ordering=genmmd", NULL};
  char* colamd[] = {"taucs.factor.ordering=colamd", NULL};
  char* amd[] = {"taucs.factor.ordering=amd", NULL};
  void* opt_arg[] = { NULL };
  
  rc = taucs_factor_solve(A,NULL,1, y,b,metis,opt_arg);
  if (rc != TAUCS_SUCCESS) return rc;
  if (rnorm(A,x,y,z)) return TAUCS_ERROR;

  rc = taucs_factor_solve(A,NULL,1, y,b,genmmd,opt_arg);
  if (rc != TAUCS_SUCCESS) return rc;
  if (rnorm(A,x,y,z)) return TAUCS_ERROR;

  rc = taucs_factor_solve(A,NULL,1, y,b,amd,opt_arg);
  if (rc != TAUCS_SUCCESS) return rc;
  if (rnorm(A,x,y,z)) return TAUCS_ERROR;

  /* colamd should fail on symmetric matrices */
  rc = taucs_factor_solve(A,NULL,1, y,b,colamd,opt_arg);
  if (rc == TAUCS_SUCCESS) return TAUCS_ERROR; 

  printf("TESING SYMMETRIC ORDERINGS SUCCEDDED\n");

  return TAUCS_SUCCESS;
}

int test_spd_factorizations(taucs_ccs_matrix* A, 
			    double* x, double* y, double* b, double* z)
{
  int rc;
  char* mf[]   = {"taucs.factor.mf=true", NULL};
  char* ll[]   = {"taucs.factor.ll=true", NULL};
  char* mfmd[] = {"taucs.factor.mf=true", "taucs.maxdepth=5", NULL};
  char* llmd[] = {"taucs.factor.ll=true", "taucs.maxdepth=5", NULL};
  char* ooc[]  = {"taucs.ooc=true", "taucs.ooc.basename=/tmp/taucs-test", NULL};
  void* opt_arg[] = { NULL };
  
  rc = taucs_factor_solve(A,NULL,1, y,b,ooc,opt_arg);
  if (rc != TAUCS_SUCCESS) return rc;
  if (rnorm(A,x,y,z)) return TAUCS_ERROR;

  rc = taucs_factor_solve(A,NULL,1, y,b,mf,opt_arg);
  if (rc != TAUCS_SUCCESS) return rc;
  if (rnorm(A,x,y,z)) return TAUCS_ERROR;

  rc = taucs_factor_solve(A,NULL,1, y,b,ll,opt_arg);
  if (rc != TAUCS_SUCCESS) return rc;
  if (rnorm(A,x,y,z)) return TAUCS_ERROR;

  /* low depth should fail  */
  rc = taucs_factor_solve(A,NULL,1, y,b,mfmd,opt_arg);
  if (rc == TAUCS_SUCCESS) return TAUCS_ERROR; 

  rc = taucs_factor_solve(A,NULL,1, y,b,llmd,opt_arg);
  if (rc == TAUCS_SUCCESS) return TAUCS_ERROR; 

  printf("TESING SPD FACTORIZATIONS SUCCEDDED\n");

  return TAUCS_SUCCESS;
}

int test_spd_factorsolve(taucs_ccs_matrix* A, 
			 double* x, double* y, double* b, double* z)
{
  int rc;
  void* F = NULL;
  char* factor[] = {"taucs.solve=false", NULL};
  char* solve [] = {"taucs.factor=false", NULL};
  void* opt_arg[] = { NULL };
  
  /* solve without a factorization should fail */
  rc = taucs_factor_solve(A,NULL,1, y,b,solve,opt_arg);
  if (rc == TAUCS_SUCCESS) return TAUCS_ERROR;

  /* solve without a factorization should fail */
  rc = taucs_factor_solve(A,&F,1, y,b,solve,opt_arg);
  if (rc == TAUCS_SUCCESS) return TAUCS_ERROR;

  rc = taucs_factor_solve(A,&F,1, y,b,factor,opt_arg);
  if (rc != TAUCS_SUCCESS) return rc;

  rc = taucs_factor_solve(A,&F,1, y,b,solve,opt_arg);
  if (rc != TAUCS_SUCCESS) return rc;
  if (rnorm(A,x,y,z)) return TAUCS_ERROR;

  printf("TESING SPD FACTORSOLVE SUCCEDDED\n");

  return TAUCS_SUCCESS;
}

int main()
{
  int xyz = 20;
  int i;

  taucs_ccs_matrix*  A;
  void* f;
  void* fc;

  double*      X;
  double*      B;
  double*      Y;
  double*      Z;

  taucs_logfile("stdout");

  A = taucs_ccs_generate_mesh3d(xyz,xyz,xyz);
  if (!A) {
    taucs_printf("Matrix generation failed\n");
    return 1;
  }

  X =(double*)malloc((A->n)*sizeof(double));
  B =(double*)malloc((A->n)*sizeof(double));
  Y =(double*)malloc((A->n)*sizeof(double));
  Z =(double*)malloc((A->n)*sizeof(double));
  if (!X || !B || !Y || !Z) {
    taucs_printf("Vector allocation failed\n");
    return 1;
  }

  for(i=0; i<A->n; i++) X[i]=((double)random()/RAND_MAX);
  taucs_ccs_times_vec(A,X,B);

  if (test_spd_factorsolve(A,X,Y,B,Z)) {
    printf("SPD FACTOR-SOLVE FAILED\n");
    return 1;
  }

  if (test_spd_orderings(A,X,Y,B,Z)) {
    printf("SPD ORDERING FAILED\n");
    return 1;
  }

  if (test_spd_factorizations(A,X,Y,B,Z)) {
    printf("SPD FACTORIZATION FAILED\n");
    return 1;
  }


  taucs_printf("test succeeded\n");
  return 0;
}
