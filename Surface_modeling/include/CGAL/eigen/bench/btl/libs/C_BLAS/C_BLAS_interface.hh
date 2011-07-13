//=====================================================
// File   :  C_BLAS_interface.hh
// Author :  L. Plagne <laurent.plagne@edf.fr)>
// Copyright (C) EDF R&D,  lun sep 30 14:23:28 CEST 2002
//=====================================================
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation; either version 2
// of the License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
//
#ifndef C_BLAS_PRODUIT_MATRICE_VECTEUR_HH
#define C_BLAS_PRODUIT_MATRICE_VECTEUR_HH

#include "f77_interface.hh"
#include <complex>
extern "C"
{
#include "cblas.h"

// #ifdef PUREBLAS
#include "blas.h"
// #endif

// void sgemm_(const char *transa, const char *transb, const int *m, const int *n, const int *k,
//            const float *alpha, const float *a, const int *lda, const float *b, const int *ldb,
//            const float *beta, float *c, const int *ldc);
//
// void sgemv_(const char *trans, const int *m, const int *n, const float *alpha,
//            const float *a, const int *lda, const float *x, const int *incx,
//            const float *beta, float *y, const int *incy);
//
// void ssymv_(const char *trans, const char* uplo,
//             const int* N, const float* alpha, const float *A,
//             const int* lda, const float *X, const int* incX,
//             const float* beta, float *Y, const int* incY);
//
// void sscal_(const int *n, const float *alpha, const float *x, const int *incx);
//
// void saxpy_(const int *n, const float *alpha, const float *x, const int *incx,
//             float *y, const int *incy);
//
// void strsv_(const char *uplo, const char *trans, const char *diag, const int *n,
//            const float *a, const int *lda, float *x, const int *incx);
//
// void scopy_(const int *n, const float *x, const int *incx, float *y, const int *incy);

  // Cholesky Factorization
// #include "mkl_lapack.h"
//   void spotrf_(const char* uplo, const int* n, float *a, const int* ld, int* info);
//   void dpotrf_(const char* uplo, const int* n, double *a, const int* ld, int* info);
  void ssytrd_(char *uplo, const int *n, float *a, const int *lda, float *d, float *e, float *tau, float *work, int *lwork, int *info );
  void sgehrd_( const int *n, int *ilo, int *ihi, float *a, const int *lda, float *tau, float *work, int *lwork, int *info );

  // LU row pivoting
//   void dgetrf_( int *m, int *n, double *a, int *lda, int *ipiv, int *info );
//   void sgetrf_(const int* m, const int* n, float *a, const int* ld, int* ipivot, int* info);
  // LU full pivoting
  void sgetc2_(const int* n, float *a, const int *lda, int *ipiv, int *jpiv, int*info );
#ifdef HAS_LAPACK
#endif
}

#define MAKE_STRING2(S) #S
#define MAKE_STRING(S) MAKE_STRING2(S)

template<class real>
class C_BLAS_interface : public f77_interface_base<real>
{
public :

  typedef typename f77_interface_base<real>::gene_matrix gene_matrix;
  typedef typename f77_interface_base<real>::gene_vector gene_vector;

  static inline std::string name( void )
  {
    return MAKE_STRING(CBLASNAME);
  }

  static  inline void matrix_vector_product(gene_matrix & A, gene_vector & B, gene_vector & X, int N)
  {
    cblas_dgemv(CblasColMajor,CblasNoTrans,N,N,1.0,A,N,B,1,0.0,X,1);
  }

  static  inline void atv_product(gene_matrix & A, gene_vector & B, gene_vector & X, int N)
  {
    cblas_dgemv(CblasColMajor,CblasTrans,N,N,1.0,A,N,B,1,0.0,X,1);
  }

  static  inline void symv(gene_matrix & A, gene_vector & B, gene_vector & X, int N)
  {
    cblas_dsymv(CblasColMajor,CblasLower,CblasTrans,N,N,1.0,A,N,B,1,0.0,X,1);
  }

  static  inline void matrix_matrix_product(gene_matrix & A, gene_matrix & B, gene_matrix & X, int N){
    cblas_dgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,N,N,N,1.0,A,N,B,N,0.0,X,N);
  }

  static  inline void transposed_matrix_matrix_product(gene_matrix & A, gene_matrix & B, gene_matrix & X, int N){
    cblas_dgemm(CblasColMajor,CblasTrans,CblasTrans,N,N,N,1.0,A,N,B,N,0.0,X,N);
  }

  static  inline void ata_product(gene_matrix & A, gene_matrix & X, int N){
    cblas_dgemm(CblasColMajor,CblasTrans,CblasNoTrans,N,N,N,1.0,A,N,A,N,0.0,X,N);
  }

  static  inline void aat_product(gene_matrix & A, gene_matrix & X, int N){
    cblas_dgemm(CblasColMajor,CblasNoTrans,CblasTrans,N,N,N,1.0,A,N,A,N,0.0,X,N);
  }

  static  inline void axpy(real coef, const gene_vector & X, gene_vector & Y, int N){
    cblas_daxpy(N,coef,X,1,Y,1);
  }

  static inline void axpby(real a, const gene_vector & X, real b, gene_vector & Y, int N){
    cblas_dscal(N,b,Y,1);
    cblas_daxpy(N,a,X,1,Y,1);
  }

};

static float fone = 1;
static float fzero = 0;
static char notrans = 'N';
static char trans = 'T';
static char nonunit = 'N';
static char lower = 'L';
static char right = 'R';
static char left = 'L';
static int intone = 1;

template<>
class C_BLAS_interface<float> : public f77_interface_base<float>
{

public :

  static inline std::string name( void )
  {
    return MAKE_STRING(CBLASNAME);
  }

  static inline void matrix_vector_product(gene_matrix & A, gene_vector & B, gene_vector & X, int N){
    #ifdef PUREBLAS
    sgemv_(&notrans,&N,&N,&fone,A,&N,B,&intone,&fzero,X,&intone);
    #else
    cblas_sgemv(CblasColMajor,CblasNoTrans,N,N,1.0,A,N,B,1,0.0,X,1);
    #endif
  }

  static inline void symv(gene_matrix & A, gene_vector & B, gene_vector & X, int N){
    #ifdef PUREBLAS
    ssymv_(&lower, &N,&fone,A,&N,B,&intone,&fzero,X,&intone);
    #else
    cblas_ssymv(CblasColMajor,CblasLower,N,1.0,A,N,B,1,0.0,X,1);
    #endif
  }

  static inline void syr2(gene_matrix & A, gene_vector & B, gene_vector & X, int N){
    #ifdef PUREBLAS
    ssyr2_(&lower,&N,&fone,B,&intone,X,&intone,A,&N);
    #else
    cblas_ssyr2(CblasColMajor,CblasLower,N,1.0,B,1,X,1,A,N);
    #endif
  }

  static inline void ger(gene_matrix & A, gene_vector & X, gene_vector & Y, int N){
    #ifdef PUREBLAS
    sger_(&N,&N,&fone,X,&intone,Y,&intone,A,&N);
    #else
    cblas_sger(CblasColMajor,N,N,1.0,X,1,Y,1,A,N);
    #endif
  }

  static inline void rot(gene_vector & A,  gene_vector & B, float c, float s, int N){
    #ifdef PUREBLAS
    srot_(&N,A,&intone,B,&intone,&c,&s);
    #else
    cblas_srot(N,A,1,B,1,c,s);
    #endif
  }

  static inline void atv_product(gene_matrix & A, gene_vector & B, gene_vector & X, int N){
    #ifdef PUREBLAS
    sgemv_(&trans,&N,&N,&fone,A,&N,B,&intone,&fzero,X,&intone);
    #else
    cblas_sgemv(CblasColMajor,CblasTrans,N,N,1.0,A,N,B,1,0.0,X,1);
    #endif
  }

  static inline void matrix_matrix_product(gene_matrix & A, gene_matrix & B, gene_matrix & X, int N){
    #ifdef PUREBLAS
    sgemm_(&notrans,&notrans,&N,&N,&N,&fone,A,&N,B,&N,&fzero,X,&N);
    #else
    cblas_sgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,N,N,N,1.0,A,N,B,N,0.0,X,N);
    #endif
  }

  static inline void transposed_matrix_matrix_product(gene_matrix & A, gene_matrix & B, gene_matrix & X, int N){
    #ifdef PUREBLAS
    sgemm_(&notrans,&notrans,&N,&N,&N,&fone,A,&N,B,&N,&fzero,X,&N);
    #else
    cblas_sgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,N,N,N,1.0,A,N,B,N,0.0,X,N);
    #endif
  }

  static inline void ata_product(gene_matrix & A, gene_matrix & X, int N){
    #ifdef PUREBLAS
    sgemm_(&trans,&notrans,&N,&N,&N,&fone,A,&N,A,&N,&fzero,X,&N);
    #else
    cblas_sgemm(CblasColMajor,CblasTrans,CblasNoTrans,N,N,N,1.0,A,N,A,N,0.0,X,N);
    #endif
  }

  static inline void aat_product(gene_matrix & A, gene_matrix & X, int N){
    #ifdef PUREBLAS
    sgemm_(&notrans,&trans,&N,&N,&N,&fone,A,&N,A,&N,&fzero,X,&N);
    #else
    cblas_sgemm(CblasColMajor,CblasNoTrans,CblasTrans,N,N,N,1.0,A,N,A,N,0.0,X,N);
    #endif
  }

  static inline void axpy(float coef, const gene_vector & X, gene_vector & Y, int N){
    #ifdef PUREBLAS
    saxpy_(&N,&coef,X,&intone,Y,&intone);
    #else
    cblas_saxpy(N,coef,X,1,Y,1);
    #endif
  }

  static inline void axpby(float a, const gene_vector & X, float b, gene_vector & Y, int N){
    #ifdef PUREBLAS
    sscal_(&N,&b,Y,&intone);
    saxpy_(&N,&a,X,&intone,Y,&intone);
    #else
    cblas_sscal(N,b,Y,1);
    cblas_saxpy(N,a,X,1,Y,1);
    #endif
  }

  static inline void cholesky(const gene_matrix & X, gene_matrix & C, int N){
    int N2 = N*N;
    scopy_(&N2, X, &intone, C, &intone);
    char uplo = 'L';
    int info = 0;
    spotrf_(&uplo, &N, C, &N, &info);
    if(info!=0) std::cerr << "spotrf_ error " << info << "\n";
  }

  static inline void partial_lu_decomp(const gene_matrix & X, gene_matrix & C, int N){
    int N2 = N*N;
    scopy_(&N2, X, &intone, C, &intone);
    char uplo = 'L';
    int info = 0;
    int * ipiv = (int*)alloca(sizeof(int)*N);
    sgetrf_(&N, &N, C, &N, ipiv, &info);
    if(info!=0) std::cerr << "sgetrf_ error " << info << "\n";
  }

  #ifdef HAS_LAPACK

  static inline void lu_decomp(const gene_matrix & X, gene_matrix & C, int N){
    int N2 = N*N;
    scopy_(&N2, X, &intone, C, &intone);
    char uplo = 'L';
    int info = 0;
    int * ipiv = (int*)alloca(sizeof(int)*N);
    int * jpiv = (int*)alloca(sizeof(int)*N);
    sgetc2_(&N, C, &N, ipiv, jpiv, &info);
  }



  static inline void hessenberg(const gene_matrix & X, gene_matrix & C, int N){
#ifdef PUREBLAS
    {
    int N2 = N*N;
    int inc = 1;
    scopy_(&N2, X, &inc, C, &inc);
    }
#else
    cblas_scopy(N*N, X, 1, C, 1);
#endif
    int info = 0;
    int ilo = 1;
    int ihi = N;
    int bsize = 64;
    int worksize = N*bsize;
    float* d = new float[N+worksize];
    sgehrd_(&N, &ilo, &ihi, C, &N, d, d+N, &worksize, &info);
    delete[] d;
  }

  static inline void tridiagonalization(const gene_matrix & X, gene_matrix & C, int N){
#ifdef PUREBLAS
    {
    int N2 = N*N;
    int inc = 1;
    scopy_(&N2, X, &inc, C, &inc);
    }
#else
    cblas_scopy(N*N, X, 1, C, 1);
#endif
    char uplo = 'U';
    int info = 0;
    int ilo = 1;
    int ihi = N;
    int bsize = 64;
    int worksize = N*bsize;
    float* d = new float[3*N+worksize];
    ssytrd_(&uplo, &N, C, &N, d, d+N, d+2*N, d+3*N, &worksize, &info);
    delete[] d;
  }
  #endif

  static inline void trisolve_lower(const gene_matrix & L, const gene_vector& B, gene_vector & X, int N){
    #ifdef PUREBLAS
    scopy_(&N, B, &intone, X, &intone);
    strsv_(&lower, &notrans, &nonunit, &N, L, &N, X, &intone);
    #else
    cblas_scopy(N, B, 1, X, 1);
    cblas_strsv(CblasColMajor, CblasLower, CblasNoTrans, CblasNonUnit, N, L, N, X, 1);
    #endif
  }

  static inline void trisolve_lower_matrix(const gene_matrix & L, const gene_matrix& B, gene_matrix & X, int N){
    #ifdef PUREBLAS
    scopy_(&N, B, &intone, X, &intone);
    strsm_(&right, &lower, &notrans, &nonunit, &N, &N, &fone, L, &N, X, &N);
    #else
    cblas_scopy(N, B, 1, X, 1);
    cblas_strsm(CblasColMajor, CblasRight, CblasLower, CblasNoTrans, CblasNonUnit, N, N, 1, L, N, X, N);
    #endif
  }

  static inline void trmm(gene_matrix & A, gene_matrix & B, gene_matrix & X, int N){
    #ifdef PUREBLAS
    strmm_(&left, &lower, &notrans,&nonunit, &N,&N,&fone,A,&N,B,&N);
    #else
    cblas_strmm(CblasColMajor, CblasLeft, CblasLower, CblasNoTrans,CblasNonUnit, N,N,1,A,N,B,N);
    #endif
  }

};


#endif



