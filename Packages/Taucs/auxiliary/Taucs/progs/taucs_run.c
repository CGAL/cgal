/*********************************************************/
/* TAUCS                                                 */
/* Author: Sivan Toledo                                  */
/*********************************************************/
/* 

TAUCS_CONFIG_DEFAULT OFF
TAUCS_CONFIG BASE
TAUCS_CONFIG DREAL
TAUCS_CONFIG SREAL
TAUCS_CONFIG DCOMPLEX
TAUCS_CONFIG SCOMPLEX
TAUCS_CONFIG FACTOR
TAUCS_CONFIG OOC_LLT
TAUCS_CONFIG METIS
TAUCS_CONFIG GENMMD
TAUCS_CONFIG COLAMD
TAUCS_CONFIG AMD
TAUCS_CONFIG GENERIC_COMPLEX
TAUCS_CONFIG MALLOC_STUBS
TAUCS_CONFIG MATRIX_GENERATORS
TAUCS_CONFIG MATRIX_IO
TAUCS_CONFIG AD_HOC_TEST
TAUCS_CONFIG_END

*/

#include <stdio.h>
#include <stdlib.h>
#include "taucs.h"

/* this only works for generic complex */
#define taucs_im(x)    ((x).i)
#define taucs_re(x)    ((x).r)

#define TRUE  1
#define FALSE 0

void rnorm(taucs_ccs_matrix* A, void* x, void* b, void* aux)
{
  double relerr;
  
  taucs_ccs_times_vec(A,x,aux);

  taucs_vec_axpby(A->n,A->flags,1.0,aux,-1.0,b,aux);

  relerr = taucs_vec_norm2(A->n,A->flags,aux) 
           / taucs_vec_norm2(A->n,A->flags,b);


  taucs_printf("relative 2-norm of the residual %.2e \n",relerr);
}

int main(int argc, char* argv[])
{
  int rc;
  int i;

  taucs_ccs_matrix*  A = NULL;
  void*      X;
  void*      B;
  void*      Y;
  void*      Z;
  void* opt_arg[] = { NULL };

  char* opt_ijv = NULL;
  char* opt_hb  = NULL;
  char* opt_log = "stdout";
  double opt_3d = -1.0;
  double opt_2d = -1.0;
  char*  opt_2d_type = "dirichlet";

  int opt_sreal    = 0;
  int opt_dreal    = 0;
  int opt_scomplex = 0;
  int opt_dcomplex = 0;
  int datatype     = TAUCS_DOUBLE;

  for (i=0; argv[i]; i++) {
    int understood = FALSE;
    
    understood |= taucs_getopt_boolean(argv[i],opt_arg,"taucs_run.sreal",&opt_sreal);
    understood |= taucs_getopt_boolean(argv[i],opt_arg,"taucs_run.dreal",&opt_dreal);
    understood |= taucs_getopt_boolean(argv[i],opt_arg,"taucs_run.scomplex",&opt_scomplex);
    understood |= taucs_getopt_boolean(argv[i],opt_arg,"taucs_run.dcomplex",&opt_dcomplex);

    understood |= taucs_getopt_string(argv[i],opt_arg,"taucs_run.ijv",&opt_ijv);
    understood |= taucs_getopt_string(argv[i],opt_arg,"taucs_run.hb", &opt_hb );
    understood |= taucs_getopt_string(argv[i],opt_arg,"taucs_run.log",&opt_log);
    understood |= taucs_getopt_double(argv[i],opt_arg,"taucs_run.mesh3d",&opt_3d);
    understood |= taucs_getopt_double(argv[i],opt_arg,"taucs_run.mesh2d",&opt_2d);
    understood |= taucs_getopt_string(argv[i],opt_arg,"taucs_run.mesh2d.type",&opt_2d_type);
    
    if (!understood) taucs_printf("taucs_run: illegal option [%s]\n",
				  argv[i]);
  }

  if (opt_sreal   ) datatype = TAUCS_SINGLE;
  if (opt_dreal   ) datatype = TAUCS_DOUBLE;
  if (opt_scomplex) datatype = TAUCS_SCOMPLEX;
  if (opt_dcomplex) datatype = TAUCS_DCOMPLEX;

  taucs_logfile(opt_log);

  if (opt_3d > 0) {
    A = taucs_ccs_generate_mesh3d((int)opt_3d,(int)opt_3d,(int)opt_3d);
    if (!A) {
      taucs_printf("Matrix generation failed\n");
      return 1;
    }
    datatype = TAUCS_DOUBLE;
  }

  if (opt_2d > 0) {
    A = taucs_ccs_generate_mesh2d((int)opt_2d,opt_2d_type);
    if (!A) {
      taucs_printf("Matrix generation failed\n");
      return 1;
    }
    datatype = TAUCS_DOUBLE;
  }


  if (opt_ijv) {
    switch (datatype) {
    case TAUCS_SINGLE:
      A = taucs_ccs_read_ijv (opt_ijv,TAUCS_SYMMETRIC | TAUCS_SINGLE); break;
      break;
    case TAUCS_DOUBLE:
      A = taucs_ccs_read_ijv (opt_ijv,TAUCS_SYMMETRIC | TAUCS_DOUBLE); break;
      break;
    case TAUCS_SCOMPLEX:
      A = taucs_ccs_read_ijv (opt_ijv,TAUCS_HERMITIAN | TAUCS_SCOMPLEX); break;
      break;
    case TAUCS_DCOMPLEX:
      A = taucs_ccs_read_ijv (opt_ijv,TAUCS_HERMITIAN | TAUCS_DCOMPLEX); break;
      break;
    default:
      taucs_printf("taucs_run: incorrect datatype\n");
      return 1;
      break;
    }      
  }

  if (opt_hb) {
    switch (datatype) {
    case TAUCS_SINGLE:
      A = taucs_ccs_read_hb (opt_hb, TAUCS_SINGLE); break;
      break;
    case TAUCS_DOUBLE:
      A = taucs_ccs_read_hb (opt_hb, TAUCS_DOUBLE); break;
      break;
    case TAUCS_SCOMPLEX:
      A = taucs_ccs_read_hb (opt_hb, TAUCS_SCOMPLEX); break;
      break;
    case TAUCS_DCOMPLEX:
      A = taucs_ccs_read_hb (opt_hb, TAUCS_DCOMPLEX); break;
      break;
    default:
      taucs_printf("taucs_run: incorrect datatype\n");
      return 1;
      break;
    }
    datatype = A->flags;
  }

  if (!A) {
    taucs_printf("taucs_run: there is no matrix!\n");
    return 1;
  }

  X = taucs_vec_create(A->n,A->flags);
  B = taucs_vec_create(A->n,A->flags);
  Y = taucs_vec_create(A->n,A->flags);
  Z = taucs_vec_create(A->n,A->flags);
  if (!X || !B || !Y || !Z) {
    taucs_printf("taucs_run: vector allocation failed\n");
    return 1;
  }

  for(i=0; i<A->n; i++) {
    if (datatype & TAUCS_SINGLE) 
      ((taucs_single*)X)[i]=(taucs_single) ((double)rand()/(double)RAND_MAX);
    if (datatype & TAUCS_DOUBLE) 
      ((taucs_double*)X)[i]=(taucs_double) ((double)rand()/(double)RAND_MAX);
    if (datatype & TAUCS_SCOMPLEX) {
      taucs_single   cre,cim;
      cre = (taucs_single) ((double)rand()/(double)RAND_MAX);
      cim = (taucs_single) ((double)rand()/(double)RAND_MAX);
      ((taucs_scomplex*)X)[i] = taucs_ccomplex_create(cre,cim);
    }
    if (datatype & TAUCS_DCOMPLEX) {
      taucs_single   zre,zim;
      zre = (taucs_double) ((double)rand()/(double)RAND_MAX);
      zim = (taucs_double) ((double)rand()/(double)RAND_MAX);
      ((taucs_dcomplex*)X)[i] = taucs_zcomplex_create(zre,zim);
    }
  }

  taucs_ccs_times_vec(A,X,B);
  
  rc = taucs_linsolve(A,NULL,1,Y,B,argv,opt_arg);

  rnorm(A,Y,B,Z);

  return 0;
}
