/*********************************************************/
/* TAUCS                                                 */
/* Author: Sivan Toledo                                  */
/*********************************************************/
/*
TAUCS_CONFIG_BEGIN
TAUCS_CONFIG_DEFAULT OFF
TAUCS_CONFIG DREAL
TAUCS_CONFIG AD_HOC_TEST
TAUCS_CONFIG_END

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include <taucs.h>

/*********************************************************/
/*                                                       */
/*********************************************************/

static double twonorm(int n, double* v)
{
  /*
  double norm;
  int i;

  for (i=0, norm=0.0; i<n; i++) norm += v[i]*v[i];

  norm = sqrt(norm);
  return norm;
  */

  double norm, ssq, scale, absvi;
  int i;

  if (n==1) return fabs(v[0]);

  scale = 0.0;
  ssq   = 1.0;

  for (i=0; i<n; i++) {
    if ( v[i] != 0 ) {
      absvi = fabs(v[i]);
      if (scale < absvi) {
	ssq   = 1.0 + ssq * (scale/absvi)*(scale/absvi);
	scale = absvi;
      } else
	ssq   = ssq + (absvi/scale)*(absvi/scale);
    }
  }
  return scale * sqrt( ssq );
}

/*********************************************************/
/* test the ordering routine                             */
/*********************************************************/

void test_ccs_order(taucs_ccs_matrix* A, char* which)
{
  int* perm;
  int* invperm;
  int  t,i;

  taucs_printf("TESTING taucs_ccs_order(...,\"%s\") FOR MEMORY LEAKS\n",
	       which);

  taucs_allocation_mark_clean();
  
  taucs_ccs_order(A,&perm,&invperm,which);
  taucs_FREE(perm); taucs_FREE(invperm);
  
  t = taucs_allocation_attempts();
  taucs_allocation_assert_clean();

  for (i=0; i<t; i++) {
    taucs_allocation_mark_clean();
    taucs_allocation_induce_failure(i);
    taucs_ccs_order(A,&perm,&invperm,which);
    if (perm || invperm) {
      taucs_printf("TAUCS ERROR REPORTING PROBLEM (%s:%d, i=%d)\n",
		   __FILE__,__LINE__,i);
      exit(1);
    }
    taucs_allocation_assert_clean();
  }

}

/*********************************************************/
/* test the factorization routine                        */
/*********************************************************/

void test_ccs_factor(taucs_ccs_matrix* A, char* which)
{
  void* L;
  int i,t;

  taucs_printf("TESTING taucs_ccs_factor(...,\"%s\") FOR MEMORY LEAKS\n",
	       which);

  taucs_allocation_mark_clean();
  
  if (!strcmp(which,"mf")) {
    L=taucs_ccs_factor_llt_mf(A);
    taucs_supernodal_factor_free(L); 
  } else if (!strcmp(which,"ll")) {
    L=taucs_ccs_factor_llt_ll(A);
    taucs_supernodal_factor_free(L); 
  } else if (!strcmp(which,"ccs ldlt")) {
    L=taucs_ccs_factor_ldlt(A);
    taucs_ccs_free(L); 
  } else if (!strcmp(which,"ccs llt")) {
    L=taucs_ccs_factor_llt(A,0.0,0);
    taucs_ccs_free(L); 
  }
  
  t = taucs_allocation_attempts();
  taucs_allocation_assert_clean();

  for (i=0; i<t; i++) {
    taucs_allocation_mark_clean();
    taucs_allocation_induce_failure(i);

    printf(">>> %d/%d\n",i,t);

    if (!strcmp(which,"mf")) {
      L=taucs_ccs_factor_llt_mf(A);
    } else if (!strcmp(which,"ll")) {
      L=taucs_ccs_factor_llt_ll(A);
    } else if (!strcmp(which,"ccs ldlt")) {
      L=taucs_ccs_factor_ldlt(A);
    } else if (!strcmp(which,"ccs llt")) {
      L=taucs_ccs_factor_llt(A,0.0,0);
    }

    if (L) {
      taucs_printf("TAUCS ERROR REPORTING PROBLEM (%s:%d, i=%d)\n",
		   __FILE__,__LINE__,i);
      exit(1);
    }
    taucs_allocation_assert_clean();
  }
}

/*********************************************************/
/* test the factorization routine                        */
/*********************************************************/

void test_ccs_solve(taucs_ccs_matrix* A, char* which, 
		    void* X, void* B)
{
  void* L;
  int i,t,e;

  taucs_printf("TESTING taucs_*_solve_*(...,\"%s\") FOR MEMORY LEAKS\n",
	       which);

  if (!strcmp(which,"sn llt")) {
    L=taucs_ccs_factor_llt_mf(A);
  } else if (!strcmp(which,"ccs llt")) {
    L=taucs_ccs_factor_llt(A,0.0,0);
  } else if (!strcmp(which,"ccs ldlt")) {
    L=taucs_ccs_factor_ldlt(A);
  }

  taucs_allocation_mark_clean();
  
  if (!strcmp(which,"sn llt")) {
    e = taucs_supernodal_solve_llt(L,X,B);
  } else if (!strcmp(which,"ccs llt")) {
    e = taucs_ccs_solve_llt(L,X,B);
  } else if (!strcmp(which,"ccs ldlt")) {
    e = taucs_ccs_solve_ldlt(L,X,B);
  }

  t = taucs_allocation_attempts();
  taucs_allocation_assert_clean();

  for (i=0; i<t; i++) {
    taucs_allocation_mark_clean();
    taucs_allocation_induce_failure(i);

    if (!strcmp(which,"sn llt")) {
      e = taucs_supernodal_solve_llt(L,X,B);
    } else if (!strcmp(which,"ccs llt")) {
      e = taucs_ccs_solve_llt(L,X,B);
    } else if (!strcmp(which,"ccs ldlt")) {
      e = taucs_ccs_solve_ldlt(L,X,B);
    }

    if (e != -1) {
      taucs_printf("TAUCS ERROR REPORTING PROBLEM (%s:%d, i=%d)\n",
		   __FILE__,__LINE__,i);
      exit(1);
    }
    taucs_allocation_assert_clean();
  }

  if (!strcmp(which,"sn llt")) {
    taucs_supernodal_factor_free(L); 
  } else if (!strcmp(which,"ccs llt")) {
    taucs_ccs_free(L); 
  } else if (!strcmp(which,"ccs ldlt")) {
    taucs_ccs_free(L); 
  }
}

/*********************************************************/
/*                                                       */
/*********************************************************/

main(int argc, char* argv[])
{
  double wtime_order;
  double wtime_permute;
  double wtime_factor;
  double wtime_solve;
  double wtime_precond_create;
  double ctime_factor;

  int i,j,t;
  double NormErr;

  taucs_ccs_matrix*  A = NULL;
  taucs_ccs_matrix*  PAPT;
  taucs_ccs_matrix*  L;

  double*      X = NULL;
  double*      B = NULL;
  double*      PX;
  double*      PB;
  double*      NX;
  double*      R;

  int*         perm;
  int*         invperm;

  taucs_logfile("stdout");

  /***********************************************************/
  /* Read matrix */
  /***********************************************************/

  A = NULL;

  if (argc < 2) {
    taucs_printf("usage: %s <matrixfile>\n",argv[0]);
    exit(1);
  }

  taucs_printf("TESTING taucs_ccs_read_hb FOR MEMORY LEAKS\n");

  taucs_allocation_mark_clean();
  A = taucs_ccs_read_hb (argv[1], TAUCS_DOUBLE);
  if (!A) {
    taucs_printf("matrix matrix file not found\n");
    exit(1);
  }
  t = taucs_allocation_attempts();
  taucs_ccs_free(A);
  taucs_allocation_assert_clean();

  for (i=0; i<t; i++) {
    taucs_allocation_mark_clean();
    taucs_allocation_induce_failure(i);
    A = taucs_ccs_read_hb (argv[1], TAUCS_DOUBLE);
    if (A) {
      taucs_printf("TAUCS ERROR REPORTING PROBLEM (%s:%d, x=%d)\n",
		   __FILE__,__LINE__,i);
      exit(1);
    }
    taucs_allocation_assert_clean();
  }

  A = taucs_ccs_read_hb (argv[1], TAUCS_DOUBLE);
  taucs_allocation_mark_clean();

  /***********************************************************/
  /* Create exact solution, compute right-hand-side          */
  /***********************************************************/

  X =(double*)malloc((A->n)*sizeof(double));
  B =(double*)malloc((A->n)*sizeof(double));
  NX=(double*)malloc((A->n)*sizeof(double));
  PX=(double*)malloc((A->n)*sizeof(double));
  PB=(double*)malloc((A->n)*sizeof(double));
  R =(double*)malloc((A->n)*sizeof(double));

  for(i=0; i<(A->n); i++)  X[i]=(double)random()/RAND_MAX;

  taucs_printf("TESTING taucs_ccs_times_vec FOR MEMORY LEAKS\n");

  taucs_allocation_mark_clean();

  taucs_ccs_times_vec(A,X,B);

  t = taucs_allocation_attempts();
  taucs_allocation_assert_clean();

  for (i=0; i<t; i++) {
    taucs_allocation_mark_clean();
    taucs_allocation_induce_failure(i);
    taucs_ccs_times_vec(A,X,B);
    taucs_allocation_assert_clean();
  }

  taucs_ccs_times_vec(A,X,B);
  taucs_allocation_mark_clean();

  /***********************************************************/
  /* order                                                   */
  /***********************************************************/

#if 0
  test_ccs_order(A,"genmmd");
  test_ccs_order(A,"mmd");
  test_ccs_order(A,"amd");
  test_ccs_order(A,"metis");
  test_ccs_order(A,"random");
  test_ccs_order(A,"identity");
#endif

  taucs_ccs_order(A,&perm,&invperm,"genmmd");
  taucs_allocation_mark_clean();

  /***********************************************************/
  /* permute                                                 */
  /***********************************************************/

  taucs_printf("TESTING taucs_ccs_permute_symmetrically FOR MEMORY LEAKS\n");
  
  taucs_allocation_mark_clean();

  wtime_permute = taucs_wtime();
  PAPT = taucs_ccs_permute_symmetrically(A,perm,invperm);
  wtime_permute = taucs_wtime() - wtime_permute;
  taucs_printf("\tPermute time  = % 10.3f seconds\n",wtime_permute);
  
  taucs_ccs_free(PAPT);
  t = taucs_allocation_attempts();
  taucs_allocation_assert_clean();

  for (i=0; i<t; i++) {
    taucs_allocation_mark_clean();
    taucs_allocation_induce_failure(i);

    PAPT = taucs_ccs_permute_symmetrically(A,perm,invperm);
    if (PAPT) {
      taucs_printf("TAUCS ERROR REPORTING PROBLEM (%s:%d, x=%d)\n",
		   __FILE__,__LINE__,i);
      exit(1);
    }
    taucs_allocation_assert_clean();
  }

  PAPT = taucs_ccs_permute_symmetrically(A,perm,invperm);
  taucs_allocation_mark_clean();
  
  /***********************************************************/
  /* factor                                                  */
  /***********************************************************/

#if 0
  test_ccs_factor(PAPT,"ll");
  test_ccs_factor(PAPT,"mf");
  test_ccs_factor(PAPT,"ccs llt");
  test_ccs_factor(PAPT,"ccs ldlt");
#endif
  
  wtime_factor = taucs_wtime();
  ctime_factor = taucs_ctime();
  
  L = taucs_ccs_factor_llt_mf(PAPT);
  
  wtime_factor = taucs_wtime() - wtime_factor;
  ctime_factor = taucs_ctime() - ctime_factor;
  taucs_printf("\tFactor time   = % 10.3f seconds  ",wtime_factor);
  taucs_printf("(%.3f cpu time)\n",ctime_factor);

  taucs_allocation_mark_clean();

  /***********************************************************/
  /* solve                                                   */
  /***********************************************************/

  taucs_dvec_permute(A->n,B,PB,perm);  /* just a loop, no need to test */

  test_ccs_solve(PAPT,"ccs llt", PX,PB);
  test_ccs_solve(PAPT,"ccs ldlt",PX,PB);
  test_ccs_solve(PAPT,"sn llt",  PX,PB);

  wtime_solve = taucs_wtime();

  taucs_supernodal_solve_llt(L,PX,PB); 

  wtime_solve = taucs_wtime() - wtime_solve;
  taucs_printf("\tSolve time    = % 10.3f seconds\n",wtime_solve);

  taucs_dvec_ipermute(A->n,PX,NX,perm); /* just a loop, no need to test */

  /***********************************************************/
  /* Compute norm of forward error                           */
  /***********************************************************/

  NormErr = 0.0;
  for(i=0; i<(A->n); i++) NormErr = max(NormErr,fabs((NX[i]-X[i])/X[i]));
  taucs_printf("main: max relative error = %.2e \n",NormErr); 

  taucs_ccs_times_vec(A,NX,R);
  for(i=0; i<(A->n); i++) R[i] -= B[i]; 
  /* now R stores the residual A*NX - B */
  taucs_printf("main: relative residual norm = %.2e \n",
	       twonorm((A->n),R)/twonorm((A->n),B));
  
  /***********************************************************/
  /* Exit                                                    */
  /***********************************************************/
  
} 

