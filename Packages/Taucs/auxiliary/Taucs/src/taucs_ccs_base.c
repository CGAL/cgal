/*********************************************************/
/* TAUCS                                                 */
/* Author: Sivan Toledo                                  */
/*********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include "taucs.h"

#ifndef TAUCS_CORE
#error "You must define TAUCS_CORE to compile this file"
#endif

/*********************************************************/
/* CCS                                                   */
/*********************************************************/

#ifdef TAUCS_CORE_GENERAL
void 
taucs_ccs_free(taucs_ccs_matrix* matrix)
{
  taucs_dccs_free(matrix);
}
#endif /*TAUCS_CORE_GENERAL*/

/* 
   Here the generic and type specific routines are different
   due to a historical accident, which forces us to set the
   flags again in the generic routine.
*/

#ifdef TAUCS_CORE_GENERAL
taucs_ccs_matrix* 
taucs_ccs_create(int m, int n, int nnz, int flags)
{
  taucs_ccs_matrix* A = NULL;

#ifdef TAUCS_DOUBLE_IN_BUILD
  if (flags & TAUCS_DOUBLE)
    A = taucs_dccs_create(m,n,nnz);
#endif

#ifdef TAUCS_SINGLE_IN_BUILD
  if (flags & TAUCS_SINGLE)	 
    A = taucs_sccs_create(m,n,nnz);
#endif

#ifdef TAUCS_DCOMPLEX_IN_BUILD
  if (flags & TAUCS_DCOMPLEX)	 
    A = taucs_zccs_create(m,n,nnz);
#endif

#ifdef TAUCS_SCOMPLEX_IN_BUILD
  if (flags & TAUCS_SCOMPLEX)	 
    A = taucs_cccs_create(m,n,nnz);
#endif
  

  if (A) {
    A->flags = flags;
    return A;
  } else {
    taucs_printf("taucs_ccs_create: no data type specifiedy\n");
    return NULL;
  }
}
#endif /*TAUCS_CORE_GENERAL*/

#ifndef TAUCS_CORE_GENERAL
taucs_ccs_matrix* 
taucs_dtl(ccs_create)(int m, int n, int nnz)
{
  taucs_ccs_matrix* matrix;

  matrix = (taucs_ccs_matrix*) taucs_malloc(sizeof(taucs_ccs_matrix));
  if (!matrix) { 
    taucs_printf("taucs_ccs_create: out of memory\n");
    return NULL; 
  }

#ifdef TAUCS_CORE_DOUBLE
  matrix->flags = TAUCS_DOUBLE;
#endif

#ifdef TAUCS_CORE_SINGLE
  matrix->flags = TAUCS_SINGLE;
#endif

#ifdef TAUCS_CORE_DCOMPLEX
  matrix->flags = TAUCS_DCOMPLEX;
#endif

#ifdef TAUCS_CORE_SINGLE
  matrix->flags = TAUCS_SINGLE;
#endif

  matrix->n = n;
  matrix->m = m;
  matrix->colptr = (int*)    taucs_malloc((n+1) * sizeof(int));
  matrix->rowind = (int*)    taucs_malloc(nnz   * sizeof(int));
  matrix->taucs_values = (taucs_datatype*) taucs_malloc(nnz * sizeof(taucs_datatype));
  if (!(matrix->colptr) || !(matrix->rowind) || !(matrix->taucs_values)) {
    taucs_printf("taucs_ccs_create: out of memory (n=%d, nnz=%d)\n",n,nnz);
    taucs_free(matrix->colptr); 
    taucs_free(matrix->rowind); 
    taucs_free(matrix->taucs_values);
    taucs_free (matrix);
    return NULL; 
  }

  return matrix;
} 

void 
taucs_dtl(ccs_free)(taucs_ccs_matrix* matrix)
{
  if (!matrix) return;

  taucs_free(matrix->rowind);
  taucs_free(matrix->colptr);
  taucs_free(matrix->taucs_values);
  taucs_free(matrix);
}

#endif /*TAUCS_CORE_GENERAL*/

/*********************************************************/
/*                                                       */
/*********************************************************/

