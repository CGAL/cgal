/*********************************************************/
/* TAUCS                                                 */
/* Author: Sivan Toledo                                  */
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

  double t1, t2;

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
    double t = 0;
    int    x = 10;
    int    y = x;
    int    i = 0;
    int    f = 0;
    int    c = 0;

    for (i=1; i<100; i++) {
      A = taucs_ccs_generate_mesh2d(x,opt_2d_type);
      if (!A) {
	taucs_printf("Matrix generation failed\n");
	return 1;
      }
      taucs_logfile("none");
      t1 = taucs_wtime();
      rc = taucs_factor_solve(A,NULL,0,NULL,NULL,argv,opt_arg);
      t2 = taucs_wtime();
      taucs_logfile(opt_log);

      t = t2-t1;
      
      taucs_printf("MINUTE BENCHMARK 2D: range [%d,%d], current test %d ==> %.2f seconds\n",
		   f,c,x,t);

      if (t > 55 && t < 65) {
	taucs_printf("MINUTE BENCHMARK 2D: %d in %.2f seconds\n",x,t);
	break;
      }
      if (t < 55) {
	if (x > f) f = x;
	if (x > c) { c = x; x = 2*x; continue; }
      }
      if (t > 65) if (x < c) c = x;
      x = (f + c) / 2;
      if (x == f) x++;
      if (x == c) x--;
      if (x == f || x == c) {
	taucs_printf("MINUTE BENCHMARK DID NOT CONVERGE\n");
      }
    }
  }

  return 0;
}
