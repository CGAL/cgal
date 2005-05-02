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
#error "This is a TAUCS core file: you must define a primitive data type"
#endif

#define RNDM ((double)random()/(double)RAND_MAX);

#ifndef max /*omer*/
#define max(x,y) ( ((x) > (y)) ? (x) : (y) )
#endif

/*********************************************************/
/*                                                       */
/*********************************************************/

#ifdef TAUCS_CORE_GENERAL

taucs_double
taucs_vec_norm2(int n, int flags, void* x)
{
  int one = 1;
#ifdef TAUCS_CONFIG_DREAL
  if (flags & TAUCS_DOUBLE)
    return (taucs_double) taucs_blas_name(dnrm2)(&n, x, &one);
#endif

#ifdef TAUCS_CONFIG_SREAL
  if (flags & TAUCS_SINGLE)
    return (taucs_double) taucs_blas_name(snrm2)(&n, x, &one);
#endif

#ifdef TAUCS_CONFIG_DCOMPLEX
  if (flags & TAUCS_DCOMPLEX)
    return (taucs_double) taucs_blas_name(dznrm2)(&n, x, &one);
#endif

#ifdef TAUCS_CONFIG_SCOMPLEX
  if (flags & TAUCS_SCOMPLEX)
    return (taucs_double) taucs_blas_name(scnrm2)(&n, x, &one);
#endif

  return taucs_get_nan();
}

void
taucs_vec_axpby(int n, int flags,
		taucs_double a, void* x,
		taucs_double b, void* y,
		void* axpby)
{
#ifdef TAUCS_CONFIG_DREAL
  if (flags & TAUCS_DOUBLE)
    taucs_dvec_axpby(n,
		     (taucs_double) a, (taucs_double*) x,
		     (taucs_double) b, (taucs_double*) y,
		     (taucs_double*) axpby);
#endif

#ifdef TAUCS_CONFIG_SREAL
  if (flags & TAUCS_SINGLE)
    taucs_svec_axpby(n,
		     (taucs_single) a, (taucs_single*) x,
		     (taucs_single) b, (taucs_single*) y,
		     (taucs_single*) axpby);
#endif

#ifdef TAUCS_CONFIG_DCOMPLEX
  if (flags & TAUCS_DCOMPLEX)
    taucs_zvec_axpby(n,
		     (taucs_double) a, (taucs_dcomplex*) x,
		     (taucs_double) b, (taucs_dcomplex*) y,
		     (taucs_dcomplex*) axpby);
#endif

#ifdef TAUCS_CONFIG_SCOMPLEX
  if (flags & TAUCS_SCOMPLEX)
    taucs_cvec_axpby(n,
		     (taucs_single) a, (taucs_scomplex*) x,
		     (taucs_single) b, (taucs_scomplex*) y,
		     (taucs_scomplex*) axpby);
#endif
}

void* taucs_vec_create(int n, int flags)
{
#ifdef TAUCS_CONFIG_DREAL
  if (flags & TAUCS_DOUBLE)
    return taucs_dvec_create(n);
#endif

#ifdef TAUCS_CONFIG_SREAL
  if (flags & TAUCS_SINGLE)
    return taucs_svec_create(n);
#endif

#ifdef TAUCS_CONFIG_DCOMPLEX
  if (flags & TAUCS_DCOMPLEX)
    return taucs_zvec_create(n);
#endif

#ifdef TAUCS_CONFIG_SCOMPLEX
  if (flags & TAUCS_SCOMPLEX)
    return taucs_cvec_create(n);
#endif

  return NULL;
}

void taucs_vec_permute(int n, int flags, void* v, void* pv, int p[])
{
#ifdef TAUCS_CONFIG_DREAL
  if (flags & TAUCS_DOUBLE)
    taucs_dvec_permute(n, (taucs_double*) v, (taucs_double*) pv, p);
#endif

#ifdef TAUCS_CONFIG_SREAL
  if (flags & TAUCS_SINGLE)
    taucs_svec_permute(n,  (taucs_single*) v, (taucs_single*) pv, p);
#endif

#ifdef TAUCS_CONFIG_DCOMPLEX
  if (flags & TAUCS_DCOMPLEX)
    taucs_zvec_permute(n,  (taucs_dcomplex*) v, (taucs_dcomplex*) pv, p);
#endif

#ifdef TAUCS_CONFIG_SCOMPLEX
  if (flags & TAUCS_SCOMPLEX)
    taucs_cvec_permute(n,  (taucs_scomplex*) v, (taucs_scomplex*) pv, p);
#endif
}

void taucs_vec_ipermute(int n, int flags, void* v, void* pv, int p[])
{
#ifdef TAUCS_CONFIG_DREAL
  if (flags & TAUCS_DOUBLE)
    taucs_dvec_ipermute(n, (taucs_double*) v, (taucs_double*) pv, p);
#endif

#ifdef TAUCS_CONFIG_SREAL
  if (flags & TAUCS_SINGLE)
    taucs_svec_ipermute(n,  (taucs_single*) v, (taucs_single*) pv, p);
#endif

#ifdef TAUCS_CONFIG_DCOMPLEX
  if (flags & TAUCS_DCOMPLEX)
    taucs_zvec_ipermute(n,  (taucs_dcomplex*) v, (taucs_dcomplex*) pv, p);
#endif

#ifdef TAUCS_CONFIG_SCOMPLEX
  if (flags & TAUCS_SCOMPLEX)
    taucs_cvec_ipermute(n,  (taucs_scomplex*) v, (taucs_scomplex*) pv, p);
#endif
}

#else
void*
taucs_dtl(vec_create)(int n)
{
  return (taucs_datatype*) taucs_malloc(n*sizeof(taucs_datatype));
} 

void
taucs_dtl(vec_axpby)(int n, 
		     taucs_real_datatype a, taucs_datatype* x,
		     taucs_real_datatype b, taucs_datatype* y,
		     taucs_datatype* axpby)
{
  int i;

  for (i=0; i<n; i++) {
#ifdef TAUCS_CORE_COMPLEX
    axpby[i] = taucs_complex_create(a * taucs_re(x[i]) + b * taucs_re(y[i]),
				    a * taucs_im(x[i]) + b * taucs_im(y[i]));
    /*
    taucs_re(axpby[i]) = a * taucs_re(x[i]) + b * taucs_re(y[i]);
    taucs_im(axpby[i]) = a * taucs_im(x[i]) + b * taucs_im(y[i]);
    */
#else
    axpby[i] = a * x[i] + b * y[i];
#endif
  }
} 

void
taucs_dtl(vec_permute)(int n, taucs_datatype v[], taucs_datatype pv[], int p[])
{
  int i;
  for (i=0; i<n; i++) pv[i] = v[p[i]];
} 

void
taucs_dtl(vec_ipermute)(int n, taucs_datatype pv[], taucs_datatype v[], int invp[]) 
{
  int i;
  for (i=0; i<n; i++) v[invp[i]] = pv[i];
} 

#endif
/*********************************************************/
/*                                                       */
/*********************************************************/

