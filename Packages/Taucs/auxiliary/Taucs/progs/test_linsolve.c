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
TAUCS_CONFIG INCOMPLETE_CHOL
TAUCS_CONFIG VAIDYA
TAUCS_CONFIG METIS
TAUCS_CONFIG GENMMD
TAUCS_CONFIG COLAMD
TAUCS_CONFIG AMD
TAUCS_CONFIG MALLOC_STUBS
TAUCS_CONFIG MATRIX_GENERATORS
TAUCS_CONFIG AD_HOC_TEST
TAUCS_CONFIG_END

*/

#include <stdio.h>
#include <stdlib.h>
#include "taucs.h"

int rnorm(taucs_ccs_matrix* A, void* x, void* b, void* aux)
{
  double relerr;
  
  taucs_ccs_times_vec(A,x,aux);

  taucs_vec_axpby(A->n,A->flags,1.0,aux,-1.0,b,aux);

  relerr = taucs_vec_norm2(A->n,A->flags,aux) 
           / taucs_vec_norm2(A->n,A->flags,b);


  taucs_printf("relative 2-norm of the residual %.2e \n",relerr);

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
  char* metis[]   = {"taucs.factor.LLT=true", "taucs.factor.ordering=metis",  NULL};
  char* genmmd[]  = {"taucs.factor.LLT=true", "taucs.factor.ordering=genmmd", NULL};
  char* amd[]     = {"taucs.factor.LLT=true", "taucs.factor.ordering=amd",    NULL};
  char* colamd[]  = {"taucs.factor.LLT=true", "taucs.factor.ordering=colamd", NULL};
  void* opt_arg[] = { NULL };
  int   test = 200;
  
  printf("TEST %d\n",test++);
  rc = taucs_linsolve(A,NULL,1, y,b,metis,opt_arg);
  if (rc != TAUCS_SUCCESS) return rc;
  if (rnorm(A,y,b,z)) return TAUCS_ERROR;

  printf("TEST %d\n",test++);
  rc = taucs_linsolve(A,NULL,1, y,b,genmmd,opt_arg);
  if (rc != TAUCS_SUCCESS) return rc;
  if (rnorm(A,y,b,z)) return TAUCS_ERROR;

  printf("TEST %d\n",test++);
  rc = taucs_linsolve(A,NULL,1, y,b,amd,opt_arg);
  if (rc != TAUCS_SUCCESS) return rc;
  if (rnorm(A,y,b,z)) return TAUCS_ERROR;

  /* colamd should fail on symmetric matrices */
  printf("TEST %d\n",test++);
  rc = taucs_linsolve(A,NULL,1, y,b,colamd,opt_arg);
  if (rc == TAUCS_SUCCESS) return TAUCS_ERROR; 

  printf("TESING SYMMETRIC ORDERINGS SUCCEDDED\n");

  return TAUCS_SUCCESS;
}

int test_spd_factorizations(taucs_ccs_matrix* A, 
			    double* x, double* y, double* b, double* z)
{
  int rc;
  char* mf[]   = {"taucs.factor.LLT=true", "taucs.factor.mf=true", NULL};
  char* ll[]   = {"taucs.factor.LLT=true", "taucs.factor.ll=true", NULL};
  char* mfmd[] = {"taucs.factor.LLT=true", "taucs.factor.mf=true", "taucs.maxdepth=5", NULL};
  char* llmd[] = {"taucs.factor.LLT=true", "taucs.factor.ll=true", "taucs.maxdepth=5", NULL};
  char* ooc[]  = {"taucs.factor.LLT=true", "taucs.ooc=true", "taucs.ooc.basename=taucs-test", NULL};
  void* opt_arg[] = { NULL };
  
  rc = taucs_linsolve(A,NULL,1, y,b,ooc,opt_arg);
  if (rc != TAUCS_SUCCESS) return rc;
  if (rnorm(A,y,b,z)) return TAUCS_ERROR;

  rc = taucs_linsolve(A,NULL,1, y,b,mf,opt_arg);
  if (rc != TAUCS_SUCCESS) return rc;
  if (rnorm(A,y,b,z)) return TAUCS_ERROR;

  rc = taucs_linsolve(A,NULL,1, y,b,ll,opt_arg);
  if (rc != TAUCS_SUCCESS) return rc;
  if (rnorm(A,y,b,z)) return TAUCS_ERROR;

  /* low depth should fail  */
  rc = taucs_linsolve(A,NULL,1, y,b,mfmd,opt_arg);
  if (rc == TAUCS_SUCCESS) return TAUCS_ERROR; 

  rc = taucs_linsolve(A,NULL,1, y,b,llmd,opt_arg);
  if (rc == TAUCS_SUCCESS) return TAUCS_ERROR; 

  printf("TESING SPD FACTORIZATIONS SUCCEDDED\n");

  return TAUCS_SUCCESS;
}

int test_spd_factorsolve(taucs_ccs_matrix* A, 
			 double* x, double* y, double* b, double* z)
{
  int rc;
  void* F = NULL;
  char* factor[] = {"taucs.factor.LLT=true", NULL};
  char* solve [] = {"taucs.factor=false", NULL};
  void* opt_arg[] = { NULL };
  int   test = 100;
  
  printf("TESING SPD FACTORSOLVE\n");

  /* solve without a factorization should fail */
  printf("TEST %d\n",test++);
  rc = taucs_linsolve(A,NULL,1, y,b,solve,opt_arg);
  if (rc == TAUCS_SUCCESS) return TAUCS_ERROR;

  /* solve without a factorization should fail */
  printf("TEST %d\n",test++);
  rc = taucs_linsolve(A,&F,1, y,b,solve,opt_arg);
  if (rc == TAUCS_SUCCESS) return TAUCS_ERROR;

  /* this should work, factor, solve, free */
  printf("TEST %d\n",test++);
  rc = taucs_linsolve(A,&F,1, y,b,factor,opt_arg);
  if (rc != TAUCS_SUCCESS) return rc;
  if (rnorm(A,y,b,z)) return TAUCS_ERROR;
  printf("TEST %d\n",test++);
  rc = taucs_linsolve(NULL,&F,0, NULL,NULL,factor,opt_arg);
  if (rc != TAUCS_SUCCESS) return rc;

  /* now first factor, then solve, then free */
  printf("TEST %d\n",test++);
  rc = taucs_linsolve(A,&F,0,NULL,NULL,factor,opt_arg);
  if (rc != TAUCS_SUCCESS) return rc;
  printf("TEST %d\n",test++);
  rc = taucs_linsolve(A,&F,1,y,b,solve,opt_arg);
  if (rc != TAUCS_SUCCESS) return rc;
  printf("TEST %d\n",test++);
  rc = taucs_linsolve(NULL,&F,0, NULL,NULL,factor,opt_arg);
  if (rc != TAUCS_SUCCESS) return rc;
  if (rnorm(A,y,b,z)) return TAUCS_ERROR;

  /* factor+solve for 4 rhs */
  printf("TEST %d\n",test++);
  rc = taucs_linsolve(A,&F,4, y,b,factor,opt_arg);
  if (rc != TAUCS_SUCCESS) return rc;
  if (rnorm(A,y,b,z)) return TAUCS_ERROR;
  if (rnorm(A,y+1*(A->n),b+1*(A->n),z+1*(A->n))) return TAUCS_ERROR;
  if (rnorm(A,y+2*(A->n),b+2*(A->n),z+2*(A->n))) return TAUCS_ERROR;
  if (rnorm(A,y+3*(A->n),b+3*(A->n),z+3*(A->n))) return TAUCS_ERROR;
  printf("TEST %d\n",test++);
  rc = taucs_linsolve(NULL,&F,0, NULL,NULL,factor,opt_arg);
  if (rc != TAUCS_SUCCESS) return rc;

  printf("TESING SPD FACTORSOLVE SUCCEDDED\n");

  return TAUCS_SUCCESS;
}

int main()
{
  int n;
  int xyz = 20;
  int i;

  taucs_ccs_matrix*  A;

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

  n = A->n;

  X =(double*)malloc(4*n*sizeof(double));
  B =(double*)malloc(4*n*sizeof(double));
  Y =(double*)malloc(4*n*sizeof(double));
  Z =(double*)malloc(4*n*sizeof(double));
  if (!X || !B || !Y || !Z) {
    taucs_printf("Vector allocation failed\n");
    return 1;
  }

  for(i=0; i<4*n; i++) X[i]=((double)rand()/RAND_MAX);
  taucs_ccs_times_vec(A,X    ,B);
  taucs_ccs_times_vec(A,X+1*n,B+1*n);
  taucs_ccs_times_vec(A,X+2*n,B+2*n);
  taucs_ccs_times_vec(A,X+3*n,B+3*n);

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
