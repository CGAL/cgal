/*********************************************************/
/* TAUCS                                                 */
/* Author: Sivan Toledo                                  */
/*********************************************************/
/* 

TAUCS_CONFIG_BEGIN
TAUCS_CONFIG_DEFAULT OFF
TAUCS_CONFIG BASE
TAUCS_CONFIG CILK
TAUCS_CONFIG DREAL
TAUCS_CONFIG LLT
TAUCS_CONFIG METIS
TAUCS_CONFIG MALLOC_STUBS
TAUCS_CONFIG MATRIX_GENERATORS
TAUCS_CONFIG AD_HOC_TEST
TAUCS_CONFIG_END

*/

#include <stdio.h>
#include "taucs.h"

#ifdef TAUCS_BLAS_UNDERSCORE
#define my_dnrm2 dnrm2_
#else
#define my_dnrm2 dnrm2
#endif

double my_dnrm2();

int main()
{
  int xyz = 30;
  int i;
  int* perm;
  int* invperm;
  
  taucs_ccs_matrix*  A;
  taucs_ccs_matrix*  PAPT;
  void*              L;

  double*      Xd;
  double*      Bd;
  double*      PXd;
  double*      PBd;
  double*      NXd;

  taucs_logfile("stdout");

  A = taucs_ccs_generate_mesh3d(xyz,xyz,xyz);
  if (!A) {
    taucs_printf("Matrix generation failed\n");
    return 1;
  }

  Xd =(double*)malloc((A->n)*sizeof(double));
  for(i=0; i<A->n; i++) (Xd)[i]=(float)((double)random()/RAND_MAX);
  taucs_ccs_times_vec(A,Xd,Bd);

  taucs_ccs_order(A,&perm,&invperm,"metis");
  if (!perm) {
    taucs_printf("Ordering failed\n");
    return 1;
  }

  PAPT = taucs_ccs_permute_symmetrically(A,perm,invperm);
  if (!PAPT) {
    taucs_printf("Permute rows and columns failed\n");
    return 1;
  }

  L = taucs_ccs_factor_llt_mf(PAPT);
  if (!L) {
    taucs_printf("Factorization failed\n");
    return 1;
  }

  taucs_vec_permute(A->n,A->flags,Bd,PBd,perm);
  taucs_supernodal_solve_llt(L,PBd,NXd); /* direct solver */
  taucs_vec_ipermute(A->n,A->flags,PXd,NXd,perm);

  {
    int one = 1;
    int i;
  
    for(i=0; i<A->n; i++) PXd[i] = NXd[i]-Xd[i];
    taucs_printf("2-norm relative error %.2e \n",
		 my_dnrm2(&(A->n),PXd,&one)/my_dnrm2(&(A->n),Xd,&one)); 
  }

  taucs_printf("test succeeded\n");
  return 0;
}
