
/*********************************************************/
/* TAUCS                                                 */
/* Author: Sivan Toledo                                  */
/*                                                       */
/* Simple complex arithmetic routines.                   */
/* They are called if the compiler does not support      */
/* complex. GCC supports complex, and so do all C99      */
/* compilers.                                            */
/*                                                       */
/*********************************************************/

#include <math.h>
#include "taucs.h"

#ifdef TAUCS_CORE_DOUBLE
double taucs_get_nan()
{
  double zero = 0.0;
  double inf  = 1.0 / zero;
  double nan  = inf - inf;
  return nan;
}
#endif

#ifdef TAUCS_CORE_DOUBLE
taucs_double taucs_dzero_const     =  0.0;
taucs_double taucs_done_const      =  1.0;
taucs_double taucs_dminusone_const = -1.0;
#endif

#ifdef TAUCS_CORE_SINGLE
taucs_single taucs_szero_const     =  0.0f;
taucs_single taucs_sone_const      =  1.0f;
taucs_single taucs_sminusone_const = -1.0f;
#endif

/*#if defined(__GNUC__) && !defined(TAUCS_CONFIG_GENERIC_COMPLEX)*/
#ifdef TAUCS_C99_COMPLEX

#ifdef TAUCS_CORE_DCOMPLEX
taucs_dcomplex taucs_zzero_const     =  0.0+0.0*_Complex_I;
taucs_dcomplex taucs_zone_const      =  1.0+0.0*_Complex_I;
taucs_dcomplex taucs_zminusone_const = -1.0+0.0*_Complex_I;
#endif

#ifdef TAUCS_CORE_SCOMPLEX
taucs_scomplex  taucs_czero_const     =  0.0f+0.0f*_Complex_I;
taucs_scomplex  taucs_cone_const      =  1.0f+0.0f*_Complex_I;
taucs_scomplex  taucs_cminusone_const = -1.0f+0.0f*_Complex_I;
#endif

#else /* TAUCS_C99_COMPLEX */

#ifdef TAUCS_CORE_DCOMPLEX
taucs_dcomplex taucs_zzero_const     = { 0.0 , 0.0 };
taucs_dcomplex taucs_zone_const      = { 1.0 , 0.0 };
taucs_dcomplex taucs_zminusone_const = {-1.0 , 0.0 };
#endif

#ifdef TAUCS_CORE_SCOMPLEX
taucs_scomplex  taucs_czero_const     = { 0.0f, 0.0f};
taucs_scomplex  taucs_cone_const      = { 1.0f, 0.0f};
taucs_scomplex  taucs_cminusone_const = {-1.0f, 0.0f};
#endif

#ifdef TAUCS_CORE_COMPLEX

taucs_datatype
taucs_dtl(complex_create_fn)(taucs_real_datatype r, taucs_real_datatype i)
{
  taucs_datatype c;
  taucs_re(c) = r;
  taucs_im(c) = i;
  return c;
}

taucs_datatype
taucs_dtl(add_fn)(taucs_datatype a, taucs_datatype b)
{
  taucs_datatype c;
  taucs_re(c) = taucs_re(a) + taucs_re(b);
  taucs_im(c) = taucs_im(a) + taucs_im(b);
  return c;
}

taucs_datatype
taucs_dtl(sub_fn)(taucs_datatype a, taucs_datatype b)
{
  taucs_datatype c;
  taucs_re(c) = taucs_re(a) - taucs_re(b);
  taucs_im(c) = taucs_im(a) - taucs_im(b);
  return c;
}

taucs_datatype
taucs_dtl(mul_fn)(taucs_datatype a, taucs_datatype b)
{
  taucs_datatype c;
  taucs_re(c) = taucs_re(a) * taucs_re(b) - taucs_im(a) * taucs_im(b);
  taucs_im(c) = taucs_re(a) * taucs_im(b) + taucs_im(a) * taucs_re(b);
  return c;
}

taucs_datatype
taucs_dtl(div_fn)(taucs_datatype a, taucs_datatype b)
{
  taucs_datatype c;
  /*double r,den; omer*/
	taucs_real_datatype r,den; 

  if (fabs(taucs_re(b)) >= fabs(taucs_im(b))) {
    r   = taucs_im(b) / taucs_re(b);
    den = taucs_re(b) + r * taucs_im(b);
    taucs_re(c) = (taucs_re(a) + r * taucs_im(a))/den;
    taucs_im(c) = (taucs_im(a) - r * taucs_re(a))/den;
  } else {
    r   = taucs_re(b) / taucs_im(b);
    den = taucs_im(b) + r * taucs_re(b);
    taucs_re(c) = (r * taucs_re(a) + taucs_im(a))/den;
    taucs_im(c) = (r * taucs_im(a) - taucs_re(a))/den;
  }
  return c;
}

taucs_datatype
taucs_dtl(conj_fn)(taucs_datatype a)
{
  taucs_datatype c;
  taucs_re(c) =   taucs_re(a);
  taucs_im(c) = - taucs_im(a);
  return c;
}

taucs_datatype
taucs_dtl(neg_fn)(taucs_datatype a)
{
  taucs_datatype c;
  taucs_re(c) = - taucs_re(a);
  taucs_im(c) = - taucs_im(a);
  return c;
}

double
taucs_dtl(abs_fn)(taucs_datatype a)
{
  double x,y,temp;

#if 1
  x = fabs(taucs_re(a));
  y = fabs(taucs_im(a));
  
  if (x==0.0) return y;
  if (y==0.0) return x;
  
  if (x > y) {
    temp = y/x;
    return ( x*sqrt(1.0+temp*temp) );
  } else {
    temp = x/y;
    return ( y*sqrt(1.0+temp*temp) );
  }
#else
  return hypot(taucs_re(a), taucs_im(a));
#endif
}

taucs_datatype
taucs_dtl(sqrt_fn)(taucs_datatype a)
{
  taucs_datatype c;
  double x,y,t;/*,w; omer*/
	taucs_real_datatype w; 

  if (taucs_re(a) == 0.0 && taucs_im(a) == 0.0) {
    taucs_re(c) = 0.0;
    taucs_im(c) = 0.0;
  } else {
    x = fabs((double) taucs_re(a));
    y = fabs((double) taucs_im(a));
    if (x >= y) {
      t = y/x;
      w = (taucs_real_datatype )(sqrt(x) * sqrt(0.5 * (1.0 + sqrt(1.0 + t * t))));
    } else {
      t = x/y;
      w = (taucs_real_datatype )(sqrt(y) * sqrt(0.5 * (t + sqrt(1.0 + t * t))));
    }

    if (taucs_re(a) > 0.0) {
      taucs_re(c) = w;
			/*taucs_im(c) = taucs_im(a) / (2.0 * w); omer*/
      taucs_im(c) = (taucs_real_datatype)(taucs_im(a) / (2.0 * w));
    } else {
      x = (taucs_im(a) >= 0.0) ? w : -w;
      taucs_im(c) = (taucs_real_datatype )x;
      /*taucs_re(c) = taucs_im(a) / (2.0 * x); omer*/
			taucs_re(c) = (taucs_real_datatype)(taucs_im(a) / (2.0 * x));
    }
  }
  
  return c;
}

#endif /* TAUCS_C99_COMPLEX */

#endif /* TAUCS_CORE_COMPLEX */






