/*********************************************************/
/* TAUCS                                                 */
/* Author: Sivan Toledo                                  */
/*********************************************************/

#include <taucs_config_tests.h>
#include <taucs_config_build.h>

/*********************************************************/
/* Cilk-related stuff                                    */
/*********************************************************/

#ifdef TAUCS_CILK
#undef TAUCS_C99_COMPLEX /* cilk2c can't process complex.h */
#endif

#ifdef TAUCS_CORE_CILK
#ifdef TAUCS_CILK
/* We are compiling a Cilk source with a Cilk compiler */


#include <cilk.h>
#include <cilk-lib.h>

#define taucs_cilk   cilk
#define taucs_spawn  spawn
#define taucs_sync   sync
#define taucs_inlet  inlet
#define taucs_Self   Self
#define taucs_Cilk_active_size Cilk_active_size

#else
/* We are compiling a Cilk source, but with a C compiler */
#define cilk
#define spawn
#define sync
#define inlet
#define Self 0
#define Cilk_active_size 1

#define taucs_cilk
#define taucs_spawn
#define taucs_sync
#define taucs_inlet
#define taucs_Self 0
#define taucs_Cilk_active_size 1
#endif
#else /* not CORE_CILK */
#define taucs_cilk
#define taucs_spawn
#define taucs_sync
#define taucs_inlet
#define taucs_Self 0
#define taucs_Cilk_active_size 1
#endif

/*********************************************************/
/* other stuff                                           */
/*********************************************************/

#ifdef TAUCS_CONFIG_DREAL
#define TAUCS_DOUBLE_IN_BUILD
#endif
#ifdef TAUCS_CONFIG_SREAL
#define TAUCS_SINGLE_IN_BUILD
#endif
#ifdef TAUCS_CONFIG_DCOMPLEX
#define TAUCS_DCOMPLEX_IN_BUILD
#endif
#ifdef TAUCS_CONFIG_SCOMPLEX
#define TAUCS_SCOMPLEX_IN_BUILD
#endif

#if   defined(TAUCS_BLAS_UNDERSCORE)
#define taucs_blas_name(x) (x##_)
#elif defined(TAUCS_BLAS_NOUNDERSCORE)
#define taucs_blas_name(x) (x)
#else
#error "taucs_blas_[no]underscore_test: linking with the BLAS failed both attempts"
#endif 

#ifdef OSTYPE_win32
typedef unsigned long ssize_t;
typedef int mode_t;
typedef int perm_t;
#define random    rand
#define srandom   srand
#endif

#define TAUCS_SUCCESS                       0
#define TAUCS_ERROR                        -1
#define TAUCS_ERROR_NOMEM                  -2
#define TAUCS_ERROR_BADARGS                -3
#define TAUCS_ERROR_INDEFINITE             -4
#define TAUCS_ERROR_MAXDEPTH               -5

#define TAUCS_INT       1024
#define TAUCS_DOUBLE    2048
#define TAUCS_SINGLE    4096
#define TAUCS_DCOMPLEX  8192
#define TAUCS_SCOMPLEX 16384

#define TAUCS_LOWER      1
#define TAUCS_UPPER      2
#define TAUCS_TRIANGULAR 4
#define TAUCS_SYMMETRIC  8
#define TAUCS_HERMITIAN  16
#define TAUCS_PATTERN    32

#define TAUCS_METHOD_LLT  1
#define TAUCS_METHOD_LDLT 2
#define TAUCS_METHOD_PLU  3

#define TAUCS_VARIANT_SNMF 1
#define TAUCS_VARIANT_SNLL 2

typedef double    taucs_double;
typedef float     taucs_single;

/* The macro TAUCS_C99_COMPLEX is defined in */
/* build/OSTYPE/taucs_config_tests if the    */
/* test program progs/taucs_c99_complex_test */
/* compiles, links, and runs.                */

/*#if defined(__GNUC__) && !defined(TAUCS_CONFIG_GENERIC_COMPLEX)*/
#ifdef TAUCS_C99_COMPLEX
/*
typedef __complex__ double taucs_dcomplex;
typedef __complex__ float  taucs_scomplex;
*/

#include <complex.h>

#undef I
#ifdef _Imaginary_I
#define TAUCS_IMAGINARY_I _Imaginary_I
#else
#define TAUCS_IMAGINARY_I _Complex_I
#endif

typedef _Complex double taucs_dcomplex;
typedef _Complex float  taucs_scomplex;

#define taucs_complex_create(r,i)  ((r)+TAUCS_IMAGINARY_I*(i))
#define taucs_ccomplex_create(r,i) ((r)+TAUCS_IMAGINARY_I*(i))
#define taucs_zcomplex_create(r,i) ((r)+TAUCS_IMAGINARY_I*(i))

#define taucs_add(x,y) ((x)+(y))
#define taucs_sub(x,y) ((x)-(y))
#define taucs_mul(x,y) ((x)*(y))
#define taucs_div(x,y) ((x)/(y))
#define taucs_neg(x)   (-(x))

#define taucs_dadd(x,y) ((x)+(y))
#define taucs_dsub(x,y) ((x)-(y))
#define taucs_dmul(x,y) ((x)*(y))
#define taucs_ddiv(x,y) ((x)/(y))
#define taucs_dneg(x)   (-(x))
#define taucs_dconj(x)  (x)
#define taucs_dimag(x)    0.0
#define taucs_dreal(x)    (x)
#define taucs_dminusone -1.0
#define taucs_done      1.0
#define taucs_dzero     0.0
#define taucs_dabs(x)   (fabs(x))
#define taucs_dsqrt(x)  (sqrt(x))

#define taucs_sadd(x,y) ((x)+(y))
#define taucs_ssub(x,y) ((x)-(y))
#define taucs_smul(x,y) ((x)*(y))
#define taucs_sdiv(x,y) ((x)/(y))
#define taucs_sneg(x)   (-(x))
#define taucs_sconj(x)  (x)
#define taucs_simag(x)    0.0f
#define taucs_sreal(x)    (x)
#define taucs_sminusone -1.0f
#define taucs_sone      1.0f
#define taucs_szero     0.0f
#define taucs_sabs(x)   ((taucs_single) fabs(x))
#define taucs_ssqrt(x)  ((taucs_single) sqrt(x))

#define taucs_zadd(x,y) ((x)+(y))
#define taucs_zsub(x,y) ((x)-(y))
#define taucs_zmul(x,y) ((x)*(y))
#define taucs_zdiv(x,y) ((x)/(y))
#define taucs_zneg(x)   (-(x))
#define taucs_zconj(x)  (conj(x))
#define taucs_zimag(x)    (cimag(x))
#define taucs_zreal(x)    (creal(x))
#define taucs_zminusone -1.0+0.0*TAUCS_IMAGINARY_I
#define taucs_zone      1.0+0.0*TAUCS_IMAGINARY_I
#define taucs_zzero     0.0+0.0*TAUCS_IMAGINARY_I
#define taucs_zabs(x)   (cabs(x))
#define taucs_zsqrt(x)  (csqrt(x))

#define taucs_cadd(x,y) ((x)+(y))
#define taucs_csub(x,y) ((x)-(y))
#define taucs_cmul(x,y) ((x)*(y))
#define taucs_cdiv(x,y) ((x)/(y))
#define taucs_cneg(x)   (-(x))
#define taucs_cconj(x)  (conjf(x))
#define taucs_cimag(x)    (cimagf(x))
#define taucs_creal(x)    (crealf(x))
#define taucs_cminusone -1.0f+0.0f*TAUCS_IMAGINARY_I
#define taucs_cone      1.0f+0.0f*TAUCS_IMAGINARY_I
#define taucs_czero     0.0f+0.0f*TAUCS_IMAGINARY_I
#define taucs_cabs(x)   (cabsf(x))
#define taucs_csqrt(x)  (csqrt(x))

#if defined(TAUCS_CORE_DOUBLE)

#define taucs_conj(x)  (x)
#define taucs_im(x)    0.0
#define taucs_re(x)    (x)
#define taucs_minusone -1.0
#define taucs_one      1.0
#define taucs_zero     0.0
#define taucs_abs(x)   (fabs(x))
#define taucs_sqrt(x)  (sqrt(x))

#elif defined(TAUCS_CORE_GENERAL)
#define taucs_im(x)    0.0
#define taucs_re(x)    (x)
#define taucs_minusone -1.0
#define taucs_one      1.0
#define taucs_zero     0.0
/*
#define taucs_conj(x)  (x)
#define taucs_abs(x)   (fabs(x))
#define taucs_sqrt(x)  (sqrt(x))
*/
#elif defined(TAUCS_CORE_SINGLE)

#define taucs_conj(x)  (x)
#define taucs_im(x)    0.0f
#define taucs_re(x)    (x)
#define taucs_minusone -1.0f
#define taucs_one      1.0f
#define taucs_zero     0.0f
#define taucs_abs(x)   ((taucs_single) fabs(x))
#define taucs_sqrt(x)  ((taucs_single) sqrt(x))

#elif defined(TAUCS_CORE_DCOMPLEX)
/*
#define taucs_conj(x)  (~(x))
#define taucs_im(x)    (__imag__ (x))
#define taucs_re(x)    (__real__ (x))
#define taucs_minusone -1.0+0.0i
#define taucs_one      1.0+0.0i
#define taucs_zero     0.0+0.0i
#define taucs_abs(x)   taucs_zabs_fn(x)
#define taucs_sqrt(x)  taucs_zsqrt_fn(x)
*/

#define taucs_conj(x)  (conj(x))
#define taucs_im(x)    (cimag(x))
#define taucs_re(x)    (creal(x))
#define taucs_minusone -1.0+0.0*TAUCS_IMAGINARY_I
#define taucs_one      1.0+0.0*TAUCS_IMAGINARY_I
#define taucs_zero     0.0+0.0*TAUCS_IMAGINARY_I
#define taucs_abs(x)   (cabs(x))
#define taucs_sqrt(x)  (csqrt(x))

#elif defined(TAUCS_CORE_SCOMPLEX)
/*
#define taucs_conj(x)  (~(x))
#define taucs_im(x)    (__imag__ (x))
#define taucs_re(x)    (__real__ (x))
#define taucs_minusone -1.0f+0.0fi
#define taucs_one      1.0f+0.0fi
#define taucs_zero     0.0f+0.0fi
#define taucs_abs(x)   taucs_cabs_fn(x)
#define taucs_sqrt(x)  taucs_csqrt_fn(x)
*/
#define taucs_conj(x)  (conjf(x))
#define taucs_im(x)    (cimagf(x))
#define taucs_re(x)    (crealf(x))
#define taucs_minusone -1.0f+0.0f*TAUCS_IMAGINARY_I
#define taucs_one      1.0f+0.0f*TAUCS_IMAGINARY_I
#define taucs_zero     0.0f+0.0f*TAUCS_IMAGINARY_I
#define taucs_abs(x)   (cabsf(x))
#define taucs_sqrt(x)  (csqrtf(x))

#endif

#else /* C99 */

typedef struct {double r,i;} taucs_dcomplex;
typedef struct {float  r,i;} taucs_scomplex;

#define taucs_zcomplex_create(r,i) taucs_zcomplex_create_fn(r,i)
#define taucs_ccomplex_create(r,i) taucs_ccomplex_create_fn(r,i)

#define taucs_dadd(x,y) ((x)+(y))
#define taucs_dsub(x,y) ((x)-(y))
#define taucs_dmul(x,y) ((x)*(y))
#define taucs_ddiv(x,y) ((x)/(y))
#define taucs_dneg(x)   (-(x))
#define taucs_dconj(x)  (x)
#define taucs_dabs(x)   (fabs(x))
#define taucs_dsqrt(x)  (sqrt(x))
#define taucs_dimag(x)   0.0
#define taucs_dreal(x)   (x)
#define taucs_dminusone -1.0
#define taucs_done     1.0
#define taucs_dzero    0.0

#define taucs_sadd(x,y) ((x)+(y))
#define taucs_ssub(x,y) ((x)-(y))
#define taucs_smul(x,y) ((x)*(y))
#define taucs_sdiv(x,y) ((x)/(y))
#define taucs_sneg(x)   (-(x))
#define taucs_sconj(x)  (x)
#define taucs_sabs(x)   ((taucs_single) fabs(x))
#define taucs_ssqrt(x)  ((taucs_single) sqrt(x))
#define taucs_sim(x)   0.0f
#define taucs_sre(x)   (x)
#define taucs_sminusone -1.0f
#define taucs_sone     1.0f
#define taucs_szero    0.0f

#define taucs_zadd(x,y) taucs_zadd_fn(x,y)
#define taucs_zsub(x,y) taucs_zsub_fn(x,y)
#define taucs_zmul(x,y) taucs_zmul_fn(x,y)
#define taucs_zdiv(x,y) taucs_zdiv_fn(x,y)
#define taucs_zneg(x)   taucs_zneg_fn(x)
#define taucs_zconj(x)  taucs_zconj_fn(x)
#define taucs_zabs(x)   taucs_zabs_fn(x)
#define taucs_zsqrt(x)  taucs_zsqrt_fn(x)
#define taucs_zimag(x)    ((x).i)
#define taucs_zreal(x)    ((x).r)
#define taucs_zminusone taucs_zminusone_const
#define taucs_zone      taucs_zone_const
#define taucs_zzero     taucs_zzero_const

#define taucs_cadd(x,y) taucs_cadd_fn(x,y)
#define taucs_csub(x,y) taucs_csub_fn(x,y)
#define taucs_cmul(x,y) taucs_cmul_fn(x,y)
#define taucs_cdiv(x,y) taucs_cdiv_fn(x,y)
#define taucs_cneg(x)   taucs_cneg_fn(x)
#define taucs_cconj(x)  taucs_cconj_fn(x)
#define taucs_cabs(x)   taucs_cabs_fn(x)
#define taucs_csqrt(x)  taucs_csqrt_fn(x)
#define taucs_cimag(x)    ((x).i)
#define taucs_creal(x)    ((x).r)
#define taucs_cminusone taucs_cminusone_const
#define taucs_cone      taucs_cone_const
#define taucs_czero     taucs_czero_const

#if defined(TAUCS_CORE_DOUBLE)

#define taucs_add(x,y) ((x)+(y))
#define taucs_sub(x,y) ((x)-(y))
#define taucs_mul(x,y) ((x)*(y))
#define taucs_div(x,y) ((x)/(y))
#define taucs_neg(x)   (-(x))
#define taucs_conj(x)  (x)
#define taucs_abs(x)   (fabs(x))
#define taucs_sqrt(x)  (sqrt(x))

#define taucs_im(x)   0.0
#define taucs_re(x)   (x)
#define taucs_minusone -1.0
#define taucs_one     1.0
#define taucs_zero    0.0

#elif defined(TAUCS_CORE_GENERAL)
/*
#define taucs_add(x,y) ((x)+(y))
#define taucs_sub(x,y) ((x)-(y))
#define taucs_mul(x,y) ((x)*(y))
#define taucs_div(x,y) ((x)/(y))
#define taucs_neg(x)   (-(x))
#define taucs_conj(x)  (x)
#define taucs_abs(x)   (fabs(x))
#define taucs_sqrt(x)  (sqrt(x))
*/
#define taucs_im(x)   0.0
#define taucs_re(x)   (x)
#define taucs_minusone -1.0
#define taucs_one     1.0
#define taucs_zero    0.0

#elif defined(TAUCS_CORE_SINGLE)

#define taucs_add(x,y) ((x)+(y))
#define taucs_sub(x,y) ((x)-(y))
#define taucs_mul(x,y) ((x)*(y))
#define taucs_div(x,y) ((x)/(y))
#define taucs_neg(x)   (-(x))
#define taucs_conj(x)  (x)
#define taucs_abs(x)   ((taucs_single) fabs(x))
#define taucs_sqrt(x)  ((taucs_single) sqrt(x))

#define taucs_im(x)   0.0f
#define taucs_re(x)   (x)
#define taucs_minusone -1.0f
#define taucs_one     1.0f
#define taucs_zero    0.0f

#elif defined(TAUCS_CORE_DCOMPLEX)

#define taucs_complex_create(r,i) taucs_zcomplex_create_fn(r,i)

#define taucs_add(x,y) taucs_zadd_fn(x,y)
#define taucs_sub(x,y) taucs_zsub_fn(x,y)
#define taucs_mul(x,y) taucs_zmul_fn(x,y)
#define taucs_div(x,y) taucs_zdiv_fn(x,y)
#define taucs_neg(x)   taucs_zneg_fn(x)
#define taucs_conj(x)  taucs_zconj_fn(x)
#define taucs_abs(x)   taucs_zabs_fn(x)
#define taucs_sqrt(x)  taucs_zsqrt_fn(x)

#define taucs_im(x)    ((x).i)
#define taucs_re(x)    ((x).r)
#define taucs_minusone taucs_zminusone_const
#define taucs_one      taucs_zone_const
#define taucs_zero     taucs_zzero_const

#elif defined(TAUCS_CORE_SCOMPLEX)

#define taucs_complex_create(r,i) taucs_ccomplex_create_fn(r,i)

#define taucs_add(x,y) taucs_cadd_fn(x,y)
#define taucs_sub(x,y) taucs_csub_fn(x,y)
#define taucs_mul(x,y) taucs_cmul_fn(x,y)
#define taucs_div(x,y) taucs_cdiv_fn(x,y)
#define taucs_neg(x)   taucs_cneg_fn(x)
#define taucs_conj(x)  taucs_cconj_fn(x)
#define taucs_abs(x)   taucs_cabs_fn(x)
#define taucs_sqrt(x)  taucs_csqrt_fn(x)

#define taucs_im(x)    ((x).i)
#define taucs_re(x)    ((x).r)
#define taucs_minusone taucs_cminusone_const
#define taucs_one      taucs_cone_const
#define taucs_zero     taucs_czero_const

#endif /* SCOMPLEX */

#endif

extern taucs_double taucs_dzero_const    ;
extern taucs_double taucs_done_const     ;
extern taucs_double taucs_dminusone_const;

extern taucs_single taucs_szero_const    ;
extern taucs_single taucs_sone_const     ;
extern taucs_single taucs_sminusone_const;

extern taucs_dcomplex taucs_zzero_const    ;
extern taucs_dcomplex taucs_zone_const     ;
extern taucs_dcomplex taucs_zminusone_const;

extern taucs_scomplex taucs_czero_const    ;
extern taucs_scomplex taucs_cone_const     ;
extern taucs_scomplex taucs_cminusone_const;

#define taucs_isnan(x) (isnan((double)(taucs_re(x))) || isnan((double)(taucs_im(x))))
#define taucs_isinf(x) (isinf((double)(taucs_re(x))) || isinf((double)(taucs_im(x))))

#ifdef TAUCS_CORE_SINGLE
#define taucs_zero_const     taucs_szero_const
#define taucs_one_const      taucs_sone_const
#define taucs_minusone_const taucs_sminusone_const

#define taucs_zero_real_const     taucs_szero_const
#define taucs_one_real_const      taucs_sone_const
#define taucs_minusone_real_const taucs_sminusone_const

#define taucs_gemm  taucs_blas_name(sgemm)
#define taucs_potrf taucs_blas_name(spotrf)
#define taucs_herk  taucs_blas_name(ssyrk)
#define taucs_trsm  taucs_blas_name(strsm)
#endif

#ifdef TAUCS_CORE_DOUBLE
#define taucs_zero_const     taucs_dzero_const
#define taucs_one_const      taucs_done_const
#define taucs_minusone_const taucs_dminusone_const

#define taucs_zero_real_const     taucs_dzero_const
#define taucs_one_real_const      taucs_done_const
#define taucs_minusone_real_const taucs_dminusone_const

#define taucs_gemm  taucs_blas_name(dgemm)
#define taucs_potrf taucs_blas_name(dpotrf)
#define taucs_herk  taucs_blas_name(dsyrk)
#define taucs_trsm  taucs_blas_name(dtrsm)
#endif

/*
#ifdef TAUCS_CORE_GENERAL
#define taucs_zero_const     taucs_dzero_const
#define taucs_one_const      taucs_done_const
#define taucs_minusone_const taucs_dminusone_const
#define taucs_zero_real_const     taucs_dzero_const
#define taucs_one_real_const      taucs_done_const
#define taucs_minusone_real_const taucs_dminusone_const
#endif
*/

#ifdef TAUCS_CORE_SCOMPLEX
#define taucs_zero_const     taucs_czero_const
#define taucs_one_const      taucs_cone_const
#define taucs_minusone_const taucs_cminusone_const

#define taucs_zero_real_const     taucs_szero_const
#define taucs_one_real_const      taucs_sone_const
#define taucs_minusone_real_const taucs_sminusone_const

#define taucs_gemm  taucs_blas_name(cgemm)
#define taucs_potrf taucs_blas_name(cpotrf)
#define taucs_herk  taucs_blas_name(cherk)
#define taucs_trsm  taucs_blas_name(ctrsm)
#endif

#ifdef TAUCS_CORE_DCOMPLEX
#define taucs_zero_const     taucs_zzero_const
#define taucs_one_const      taucs_zone_const
#define taucs_minusone_const taucs_zminusone_const

#define taucs_zero_real_const     taucs_dzero_const
#define taucs_one_real_const      taucs_done_const
#define taucs_minusone_real_const taucs_dminusone_const

#define taucs_gemm  taucs_blas_name(zgemm)
#define taucs_potrf taucs_blas_name(zpotrf)
#define taucs_herk  taucs_blas_name(zherk)
#define taucs_trsm  taucs_blas_name(ztrsm)
#endif

/*********************************************************/
/*                                                       */
/*********************************************************/

typedef struct
{
  int     n;    /* columns                      */
  int     m;    /* rows; don't use if symmetric   */
  int     flags;
  int*    colptr; /* pointers to where columns begin in rowind and values. */
                  /* 0-based. Length is (n+1). */
  int*    rowind; /* row indices */

  union {
    void*           v;
    taucs_double*   d;
    taucs_single*   s;
    taucs_dcomplex* z;
    taucs_scomplex* c;
  } values;

} taucs_ccs_matrix;

typedef struct {
  int   type;
  int   nmatrices;
  void* type_specific;

  /* the following may change! do not rely on them. */
  double nreads, nwrites, bytes_read, bytes_written, read_time, write_time;
} taucs_io_handle;

/* generate all the prototypes */

#define taucs_datatype taucs_double
#define taucs_real_datatype taucs_double
#define taucs_dtl(X) taucs_d##X
#include "taucs_private.h"
#undef taucs_real_datatype
#undef taucs_datatype
#undef taucs_dtl

#define taucs_datatype taucs_single
#define taucs_real_datatype taucs_single
#define taucs_dtl(X) taucs_s##X
#include "taucs_private.h"
#undef taucs_real_datatype
#undef taucs_datatype
#undef taucs_dtl

#define taucs_datatype taucs_dcomplex
#define taucs_real_datatype taucs_double
#define taucs_dtl(X) taucs_z##X
#include "taucs_private.h"
#undef taucs_real_datatype
#undef taucs_datatype
#undef taucs_dtl

#define taucs_datatype taucs_scomplex
#define taucs_real_datatype taucs_single
#define taucs_dtl(X) taucs_c##X
#include "taucs_private.h"
#undef taucs_real_datatype
#undef taucs_datatype
#undef taucs_dtl

/*********************************************************/
/*                                                       */
/*********************************************************/

/* now define the data type for the file that we compile now */

#ifdef TAUCS_CORE_DOUBLE
#define TAUCS_CORE
#define TAUCS_CORE_REAL
#define TAUCS_CORE_DATATYPE TAUCS_DOUBLE
typedef taucs_double taucs_datatype;
#define taucs_dtl(X) taucs_d##X
#define taucs_values values.d
#define taucs_iszero(x) ((x) == 0.0)
typedef double taucs_real_datatype; /* omer: this is the datatype of the real and imaginary part of the datatype*/
#endif

#ifdef TAUCS_CORE_GENERAL
#define TAUCS_CORE
#define TAUCS_CORE_DATATYPE TAUCS_DOUBLE
typedef taucs_double taucs_datatype;
typedef double taucs_real_datatype; 
/*
#define TAUCS_CORE_REAL
#define TAUCS_CORE_DATATYPE TAUCS_DOUBLE
#define taucs_values values.d
#define taucs_dtl(X) taucs_g##X
#define taucs_iszero(x) ((x) == 0.0)
*/
#endif

#ifdef  TAUCS_CORE_SINGLE
#define TAUCS_CORE
#define TAUCS_CORE_REAL
#define TAUCS_CORE_DATATYPE TAUCS_SINGLE
typedef taucs_single taucs_datatype;
#define taucs_dtl(X) taucs_s##X
#define taucs_values values.s
#define taucs_iszero(x) ((x) == 0.0f)
typedef float taucs_real_datatype; /* omer: this is the datatype of the real and imaginary part of the datatype*/
#endif

#ifdef  TAUCS_CORE_DCOMPLEX
#define TAUCS_CORE
#define TAUCS_CORE_COMPLEX
#define TAUCS_CORE_DATATYPE TAUCS_DCOMPLEX
typedef taucs_dcomplex taucs_datatype;
#define taucs_dtl(X) taucs_z##X
#define taucs_values values.z
#define taucs_iszero(x) (taucs_re(x) == 0.0 && taucs_im(x) == 0.0)
typedef double taucs_real_datatype; /* omer: this is the datatype of the real and imaginary part of the datatype*/
#endif

#ifdef  TAUCS_CORE_SCOMPLEX
#define TAUCS_CORE
#define TAUCS_CORE_COMPLEX
#define TAUCS_CORE_DATATYPE TAUCS_SCOMPLEX
typedef taucs_scomplex taucs_datatype;
#define taucs_dtl(X) taucs_c##X
#define taucs_values values.c
#define taucs_iszero(x) (taucs_re(x) == 0.0f && taucs_im(x) == 0.0f)
typedef float taucs_real_datatype; /* omer: this is the datatype of the real and imaginary part of the datatype*/
#endif

#ifndef TAUCS_CORE_DATATYPE
typedef taucs_double taucs_datatype;
typedef double taucs_real_datatype; 
#endif

/*********************************************************/
/*                                                       */
/*********************************************************/

double taucs_get_nan(void);

/* 
   routines for testing memory allocation.
   Mostly useful for testing programs
   that hunt for memory leaks.
*/

double taucs_allocation_amount(void);
int    taucs_allocation_count(void);
int    taucs_allocation_attempts(void);
void   taucs_allocation_assert_clean(void);
void   taucs_allocation_mark_clean(void);
void   taucs_allocation_induce_failure(int i);

/* 
   these are meant to allow allocation 
   and more importantly, deallocation,
   within the testing programs.
*/

#include <stdlib.h>

void* taucs_malloc (size_t size)              ;
void* taucs_calloc (size_t nmemb, size_t size);
void* taucs_realloc(void* ptr, size_t size)   ;
void  taucs_free   (void* ptr)                ;

#if defined(TAUCS_CORE) 

#if defined(TAUCS_MEMORY_TEST_yes)

#include <stdlib.h>

void* taucs_internal_calloc(size_t nmemb, size_t size,char* file, int line);
void* taucs_internal_malloc(size_t size,              char* file, int line);
void* taucs_internal_realloc(void *ptr, size_t size,   char* file, int line);
void  taucs_internal_free(void *ptr,                   char* file, int line);

/*
 #define realloc(x,y) taucs_internal_realloc(x,y,__FILE__,__LINE__)
 #define malloc(x)    taucs_internal_malloc(x,__FILE__,__LINE__)
 #define calloc(x,y)  taucs_internal_calloc(x,y,__FILE__,__LINE__)
 #define free(x)      taucs_internal_free(x,__FILE__,__LINE__)
*/

#define taucs_realloc(x,y) taucs_internal_realloc(x,y,__FILE__,__LINE__)
#define taucs_malloc(x)    taucs_internal_malloc(x,__FILE__,__LINE__)
#define taucs_calloc(x,y)  taucs_internal_calloc(x,y,__FILE__,__LINE__)
#define taucs_free(x)      taucs_internal_free(x,__FILE__,__LINE__)

#define realloc(x,y) taucs_must_not_call_realloc_directly(x,y)
#define malloc(x)    taucs_must_not_call_malloc_directly(x)
#define calloc(x,y)  taucs_must_not_call_calloc_directly(x,y)
#define free(x)      taucs_must_not_call_free_directly(x)

#else /* TAUCS_CORE, but not memory testing */

void* taucs_calloc_stub(size_t nmemb, size_t size);
void* taucs_malloc_stub(size_t size);
void* taucs_realloc_stub(void *ptr, size_t size);
void  taucs_free_stub(void *ptr);

#define realloc(x,y) taucs_must_not_call_realloc_directly(x,y)
#define malloc(x)    taucs_must_not_call_malloc_directly(x)
#define calloc(x,y)  taucs_must_not_call_calloc_directly(x,y)
#define free(x)      taucs_must_not_call_free_directly(x)

#define taucs_realloc(x,y) taucs_realloc_stub(x,y)
#define taucs_malloc(x)    taucs_malloc_stub(x)
#define taucs_calloc(x,y)  taucs_calloc_stub(x,y)
#define taucs_free(x)      taucs_free_stub(x)

#endif
#endif

/*********************************************************/
/*                                                       */
/*********************************************************/

#ifndef max
#define max(x,y) ( ((x) > (y)) ? (x) : (y) )
#endif

#ifndef min
#define min(x,y) ( ((x) < (y)) ? (x) : (y) )
#endif

/*********************************************************/
/*                                                       */
/*********************************************************/

/* externs */
extern int ireadhb_(char*, char*, int*, int*, int*);
extern int creadhb_(char*, int*, int*, int*, int*, int*, taucs_scomplex*);
extern int dreadhb_(char*, int*, int*, int*, int*, int*, taucs_double*);
extern int sreadhb_(char*, int*, int*, int*, int*, int*, taucs_single*);
extern int zreadhb_(char*, int*, int*, int*, int*, int*, taucs_dcomplex*);

extern int amdexa_(int*, int*, int*, int*, int*, int*, int*, int*, int*, 
			int*, int*, int*, int*, int*, int*);
extern int amdtru_(int*, int*, int*, int*, int*, int*, int*, int*, int*, 
			int*, int*, int*, int*, int*, int*);
extern int amdbar_(int*, int*, int*, int*, int*, int*, int*, int*, int*, 
			int*, int*, int*, int*, int*, int*);
extern int genmmd_(int*, int*, int*, int*, int*, int*, int*, int*, int*, 
			int*, int*, int*);

/*********************************************************/
/*                                                       */
/*********************************************************/

#if (defined(OSTYPE_irix) || defined(OSTYPE_solaris))

#include <math.h>
#include <ieeefp.h>
#define isinf(x) (!finite((x)) && !isnan((x)))

#elif defined(OSTYPE_win32)

#include <float.h>
#define isnan(x)  (_isnan(x))
#define isinf(x)  (!(_finite(x)) && !(_isnan(x)))
#define finite(x) (_finite(x))

#endif

/* If these are mactors (e.g., gcc -std=c99), do not declare   */
/* otherwise, declare them, since they are not always declared */
/* in math.h (e.g., gcc -std=c89 -pedantic); these are for     */
/* gcc 3.3.1                                                   */

#ifndef isnan
extern int isnan(double);
#endif
#ifndef finite
extern int finite(double);
#endif
#ifndef isinf
extern int isinf(double);
#endif

extern int taucs_potrf(char*, int*, taucs_datatype*, int*, int*);
extern int taucs_trsm(char *, char *, char *, char *, 
			int*, int*, taucs_datatype*, taucs_datatype*, int *, 
			taucs_datatype*, int *);
extern int taucs_gemm(char *, char *, int*, int*, int *,
			taucs_datatype*, taucs_datatype*, int *, taucs_datatype*, int *, 
			taucs_datatype*, taucs_datatype*, int*);
extern int taucs_herk(char *, char *, 
		      int *, int *, 
		      taucs_real_datatype*, 
		      taucs_datatype*, int *, 
		      taucs_real_datatype*, 
		      taucs_datatype*, int *);

taucs_double taucs_blas_name(dnrm2)(int*, taucs_double*, int*);
taucs_single taucs_blas_name(snrm2)(int*, taucs_single*, int*);
taucs_double taucs_blas_name(dznrm2)(int*, taucs_dcomplex*, int*);
taucs_single taucs_blas_name(scnrm2)(int*, taucs_scomplex*, int*);

/*********************************************************/
/*                                                       */
/*********************************************************/


