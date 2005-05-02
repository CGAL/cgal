/*********************************************************/
/* TAUCS                                                 */
/* Author: Sivan Toledo                                  */
/*********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include "taucs.h"

/*** superlu includes ***/

#include "dsp_defs.h"
#include "util.h"

typedef struct {
  SuperMatrix L,U;
  int*        perm_c;
  int*        perm_r;
} taucs_superlu_factor;

/*********************************************************/
/* calling superlu                                       */
/*********************************************************/

void* 
taucs_ccs_factor_superlu(taucs_ccs_matrix* A)
{
  double* a;
  int*    asub;
  int*    xa;
  SuperMatrix sA,sAC;
  double* b;
  int     n,m,nnz,info;
  int  i,j,ip;
  int* etree;
  taucs_superlu_factor* F;

  double wtime_prepare;
  double wtime_factor;

  wtime_prepare = taucs_wtime();

  assert(A->flags & TAUCS_SYMMETRIC);
  assert(A->flags & TAUCS_LOWER);

  n = m = A->n;
  nnz = 2 * (A->colptr)[n] - n;

  a    = (double*) taucs_malloc(nnz * sizeof(double));
  asub = (int*)    taucs_malloc(nnz * sizeof(int));
  xa   = (int*)    taucs_malloc((n+1)*sizeof(int));
  etree  = (int*) taucs_malloc(n*sizeof(int));
 
  F = (taucs_superlu_factor*) taucs_malloc(sizeof(taucs_superlu_factor));

  (F->perm_r) = (int*) taucs_malloc(n*sizeof(int));
  (F->perm_c) = (int*) taucs_malloc(n*sizeof(int));

  for (i=0; i<n; i++) (F->perm_c)[i] = 0;

  for (j=0; j<n; j++) {
    for (ip = (A->colptr)[j]; ip<(A->colptr)[j+1]; ip++) {
      i = (A->rowind)[ip];
      if (i==j) 
	(F->perm_c)[j]++;
      else {
	(F->perm_c)[i]++;
	(F->perm_c)[j]++;
      }
    }
  }
  
  xa[0] = 0;
  for (j=1; j<=n; j++) xa[j] = xa[j-1] + (F->perm_c)[j-1];
  for (j=0; j< n; j++) (F->perm_c)[j] = xa[j];

  assert(nnz = xa[n]);

  for (j=0; j<n; j++) {
    for (ip = (A->colptr)[j]; ip<(A->colptr)[j+1]; ip++) {
      i = (A->rowind)[ip];
      if (i==j) {
	asub[(F->perm_c)[j]] = i;
	a   [(F->perm_c)[j]] = (A->values)[ip];
	(F->perm_c)[j]++;
      } else {
	asub[(F->perm_c)[j]] = i;
	a   [(F->perm_c)[j]] = (A->values)[ip];
	(F->perm_c)[j]++;

	asub[(F->perm_c)[i]] = j;
	a   [(F->perm_c)[i]] = (A->values)[ip];
	(F->perm_c)[i]++;
      }
    }
  }


  for (j=0; j<n; j++) {
    for (ip=xa[j]+1; ip<xa[j+1]; ip++) {
      int    key  = asub[ip];
      double keyv = a   [ip];
      int kp = ip-1;
      while (kp >= xa[j] && asub[kp] > key) {
	asub[kp+1] = asub[kp];
	a   [kp+1] = a   [kp];
	kp--;
      }
      asub[kp+1] = key;
      a   [kp+1] = keyv;
    }
  }

#if 0
  printf("xa=[ ");
  for (j=0; j<=n; j++) printf("%d ",xa[j]);
  printf("]\n");

  printf("asub=[ ");
  for (j=0; j<nnz; j++) printf("%d ",asub[j]);
  printf("]\n");

  printf("a=[ ");
  for (j=0; j<nnz; j++) printf("%lf ",a[j]);
  printf("]\n");
#endif

  for (j=0; j< n; j++) assert((F->perm_c)[j] = xa[j+1]);

  wtime_prepare = taucs_wtime() - wtime_prepare;
  taucs_printf("\t\tSuperLU preparation time = % 10.3f seconds\n",wtime_prepare);

  wtime_factor = taucs_wtime();

  dCreate_CompCol_Matrix(&sA,m,n,nnz,a,asub,xa,NC,_D,GE);

  for (i=0; i<n; i++) (F->perm_r)[i] = (F->perm_c)[i] = i;

#if 0
  printf("perm_c=[ ");
  for (j=0; j<n; j++) printf("%d ",(F->perm_c)[j]);
  printf("]\n");
#endif

  //get_perm_c(2,&sA,(F->perm_c));
  sp_preorder("N", &sA, (F->perm_c), etree, &sAC);

#if 0
  printf("perm_c=[ ");
  for (j=0; j<n; j++) printf("%d ",(F->perm_c)[j]);
  printf("]\n");
#endif

  for (i=0; i<n; i++) (F->perm_r)[i] = (F->perm_c)[i];

#if 0
  b = (double*) taucs_malloc(n*sizeof(int));
  for (i=0; i<n; i++) b[i] = 1.0;
  dCreate_Dense_Matrix(&sB,n,1,b,n,DN,_D,GE);
  dgssv(&sA,(F->perm_c),(F->perm_r),&sL,&sU,&sB,&info);
#endif

  StatInit(100,8); /* panel size, relax */

  dgstrf("N", /* not a refactorization */
	 &sAC,
	 0.0, /* don't pivot! */
	 0.0, /* drop tolerance, not implemented */
	 8,   /* relax */
	 100, /* panel size */
	 etree,
	 NULL, /* work */
	 0,    /* lwork, no memory preallocated */
	 (F->perm_r), 
	 (F->perm_c),
	 &(F->L),
	 &(F->U),
	 &info);

  taucs_printf("\t\ttaucs_ccs_factor_superlu info=%d\n",info);

  {
    mem_usage_t usage;
    extern SuperLUStat_t SuperLUStat;

    dQuerySpace(&(F->L), &(F->U), 100,&usage);
    taucs_printf("\t\ttaucs_ccs_factor_superlu %.2e bytes for L+U\n",usage.for_lu);
    taucs_printf("\t\ttaucs_ccs_factor_superlu %.2e bytes total\n",usage.total_needed);
    taucs_printf("\t\ttaucs_ccs_factor_superlu %.2e expansions\n",(float) usage.expansions);

    taucs_printf("\t\ttaucs_ccs_factor_superlu %.2e flops\n",(float) (SuperLUStat.ops)[FACT]);
  }

  StatFree();

#if 0
  printf("perm_r=[ ");
  for (j=0; j<n; j++) printf("%d ",(F->perm_r)[j]);
  printf("]\n");
#endif

  wtime_factor = taucs_wtime() - wtime_factor;
  taucs_printf("\t\tSuperLU factor time = % 10.3f seconds\n",wtime_factor);

  return F;
}

int 
taucs_ccs_solve_superlu(void* vF, double* x, double* b)
{
  taucs_superlu_factor* F = (taucs_superlu_factor*) vF;
  int n,i;
  SuperMatrix X,B;
  int info;
  double* t;

  n = (F->L).ncol;

  t = (double*) taucs_malloc(n*sizeof(double));
  for (i=0; i<n; i++) t[i] = b[i];
  dCreate_Dense_Matrix(&B,n,1,t,n,DN,_D,GE);

  dgstrs("T", /* no transpose */
	 &(F->L),&(F->U),F->perm_r,F->perm_c,
	 &B,
	 &info
	 );

  for (i=0; i<n; i++) x[i] = t[i];

  taucs_printf("\t\ttaucs_ccs_solve_superlu info=%d\n",info);
  
  return info;
}
