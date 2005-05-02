/*********************************************************/
/* TAUCS                                                 */
/* Author: Sivan Toledo                                  */
/*********************************************************/
/*
TAUCS_CONFIG_BEGIN
TAUCS_CONFIG_DEFAULT OFF
TAUCS_CONFIG BASE
TAUCS_CONFIG DREAL
TAUCS_CONFIG DCOMPLEX
TAUCS_CONFIG COLAMD
TAUCS_CONFIG ORDERING
TAUCS_CONFIG MATRIX_IO
TAUCS_CONFIG MALLOC_STUBS
TAUCS_CONFIG OOC_LU
TAUCS_CONFIG AD_HOC_TEST
TAUCS_CONFIG_END
*/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include <taucs.h>

#ifdef TAUCS_BLAS_UNDERSCORE
#define my_dnrm2 dnrm2_
#else
#define my_dnrm2 dnrm2
#endif

/*********************************************************/
/*                                                       */
/*********************************************************/

int 
main(int argc, char* argv[])
{
  void* x;
  void* b;
  taucs_ccs_matrix* A;
  void* LU;
  char  fname[256];
  double memory_mb = -1.0;
  int    mb = -1;
  int* perm;
  int* invperm;

  taucs_logfile("stdout");

  if (argc < 2) {
    fprintf(stderr,"usage: %s filename\n",argv[0]);
    fprintf(stderr,"       filename-A.bin is the matrix\n");
    fprintf(stderr,"       filename-b.bin is the rhs\n");
    fprintf(stderr,"       filename-x.bin will contain the solution\n");
    exit(1);
  }

  sprintf(fname,"%s-A.bin",argv[1]);
  A = taucs_ccs_read_binary(fname);

  if (argc > 2) {
    sscanf(argv[2],"%d",&mb);
  }

  sprintf(fname,"%s-b.bin",argv[1]);
  b = taucs_vec_read_binary(A->m, A->flags, fname);

  taucs_ccs_order(A,&perm,&invperm,"colamd");

  if (mb > 0)
    memory_mb = (double) mb;
  else
    memory_mb = ((double) (-mb)) * taucs_available_memory_size()/1048576.0;

  LU = taucs_io_create_multifile(argv[1]);
  taucs_ooc_factor_lu(A, perm, LU, memory_mb*1048576.0);

  if (A->flags & TAUCS_DOUBLE) 
    x = malloc((A->n) * sizeof(taucs_double));
  else if (A->flags & TAUCS_DCOMPLEX) 
    x = malloc((A->n) * sizeof(taucs_dcomplex));
  else
    assert(0);

  taucs_ooc_solve_lu(LU, x, b);

  sprintf(fname,"%s-x.bin",argv[1]);
  taucs_vec_write_binary(A->n, A->flags, x, fname);

  free(perm);
  free(invperm);

  return 0;
}
