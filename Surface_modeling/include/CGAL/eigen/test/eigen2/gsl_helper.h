// This file is part of Eigen, a lightweight C++ template library
// for linear algebra. Eigen itself is part of the KDE project.
//
// Copyright (C) 2008 Gael Guennebaud <g.gael@free.fr>
//
// Eigen is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 3 of the License, or (at your option) any later version.
//
// Alternatively, you can redistribute it and/or
// modify it under the terms of the GNU General Public License as
// published by the Free Software Foundation; either version 2 of
// the License, or (at your option) any later version.
//
// Eigen is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License or the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License and a copy of the GNU General Public License along with
// Eigen. If not, see <http://www.gnu.org/licenses/>.

#ifndef EIGEN_GSL_HELPER
#define EIGEN_GSL_HELPER

#include <Eigen/Core>

#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>

namespace Eigen {

template<typename Scalar, bool IsComplex = NumTraits<Scalar>::IsComplex> struct GslTraits
{
  typedef gsl_matrix* Matrix;
  typedef gsl_vector* Vector;
  static Matrix createMatrix(int rows, int cols) { return gsl_matrix_alloc(rows,cols); }
  static Vector createVector(int size) { return gsl_vector_alloc(size); }
  static void free(Matrix& m) { gsl_matrix_free(m); m=0; }
  static void free(Vector& m) { gsl_vector_free(m); m=0; }
  static void prod(const Matrix& m, const Vector& v, Vector& x) { gsl_blas_dgemv(CblasNoTrans,1,m,v,0,x); }
  static void cholesky(Matrix& m) { gsl_linalg_cholesky_decomp(m); }
  static void cholesky_solve(const Matrix& m, const Vector& b, Vector& x) { gsl_linalg_cholesky_solve(m,b,x); }
  static void eigen_symm(const Matrix& m, Vector& eval, Matrix& evec)
  {
    gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc(m->size1);
    Matrix a = createMatrix(m->size1, m->size2);
    gsl_matrix_memcpy(a, m);
    gsl_eigen_symmv(a,eval,evec,w);
    gsl_eigen_symmv_sort(eval, evec, GSL_EIGEN_SORT_VAL_ASC);
    gsl_eigen_symmv_free(w);
    free(a);
  }
  static void eigen_symm_gen(const Matrix& m, const Matrix& _b, Vector& eval, Matrix& evec)
  {
    gsl_eigen_gensymmv_workspace * w = gsl_eigen_gensymmv_alloc(m->size1);
    Matrix a = createMatrix(m->size1, m->size2);
    Matrix b = createMatrix(_b->size1, _b->size2);
    gsl_matrix_memcpy(a, m);
    gsl_matrix_memcpy(b, _b);
    gsl_eigen_gensymmv(a,b,eval,evec,w);
    gsl_eigen_symmv_sort(eval, evec, GSL_EIGEN_SORT_VAL_ASC);
    gsl_eigen_gensymmv_free(w);
    free(a);
  }
};

template<typename Scalar> struct GslTraits<Scalar,true>
{
  typedef gsl_matrix_complex* Matrix;
  typedef gsl_vector_complex* Vector;
  static Matrix createMatrix(int rows, int cols) { return gsl_matrix_complex_alloc(rows,cols); }
  static Vector createVector(int size) { return gsl_vector_complex_alloc(size); }
  static void free(Matrix& m) { gsl_matrix_complex_free(m); m=0; }
  static void free(Vector& m) { gsl_vector_complex_free(m); m=0; }
  static void cholesky(Matrix& m) { gsl_linalg_complex_cholesky_decomp(m); }
  static void cholesky_solve(const Matrix& m, const Vector& b, Vector& x) { gsl_linalg_complex_cholesky_solve(m,b,x); }
  static void prod(const Matrix& m, const Vector& v, Vector& x)
  { gsl_blas_zgemv(CblasNoTrans,gsl_complex_rect(1,0),m,v,gsl_complex_rect(0,0),x); }
  static void eigen_symm(const Matrix& m, gsl_vector* &eval, Matrix& evec)
  {
    gsl_eigen_hermv_workspace * w = gsl_eigen_hermv_alloc(m->size1);
    Matrix a = createMatrix(m->size1, m->size2);
    gsl_matrix_complex_memcpy(a, m);
    gsl_eigen_hermv(a,eval,evec,w);
    gsl_eigen_hermv_sort(eval, evec, GSL_EIGEN_SORT_VAL_ASC);
    gsl_eigen_hermv_free(w);
    free(a);
  }
  static void eigen_symm_gen(const Matrix& m, const Matrix& _b, gsl_vector* &eval, Matrix& evec)
  {
    gsl_eigen_genhermv_workspace * w = gsl_eigen_genhermv_alloc(m->size1);
    Matrix a = createMatrix(m->size1, m->size2);
    Matrix b = createMatrix(_b->size1, _b->size2);
    gsl_matrix_complex_memcpy(a, m);
    gsl_matrix_complex_memcpy(b, _b);
    gsl_eigen_genhermv(a,b,eval,evec,w);
    gsl_eigen_hermv_sort(eval, evec, GSL_EIGEN_SORT_VAL_ASC);
    gsl_eigen_genhermv_free(w);
    free(a);
  }
};

template<typename MatrixType>
void convert(const MatrixType& m, gsl_matrix* &res)
{
//   if (res)
//     gsl_matrix_free(res);
  res = gsl_matrix_alloc(m.rows(), m.cols());
  for (int i=0 ; i<m.rows() ; ++i)
    for (int j=0 ; j<m.cols(); ++j)
      gsl_matrix_set(res, i, j, m(i,j));
}

template<typename MatrixType>
void convert(const gsl_matrix* m, MatrixType& res)
{
  res.resize(int(m->size1), int(m->size2));
  for (int i=0 ; i<res.rows() ; ++i)
    for (int j=0 ; j<res.cols(); ++j)
      res(i,j) = gsl_matrix_get(m,i,j);
}

template<typename VectorType>
void convert(const VectorType& m, gsl_vector* &res)
{
  if (res) gsl_vector_free(res);
  res = gsl_vector_alloc(m.size());
  for (int i=0 ; i<m.size() ; ++i)
      gsl_vector_set(res, i, m[i]);
}

template<typename VectorType>
void convert(const gsl_vector* m, VectorType& res)
{
  res.resize (m->size);
  for (int i=0 ; i<res.rows() ; ++i)
    res[i] = gsl_vector_get(m, i);
}

template<typename MatrixType>
void convert(const MatrixType& m, gsl_matrix_complex* &res)
{
  res = gsl_matrix_complex_alloc(m.rows(), m.cols());
  for (int i=0 ; i<m.rows() ; ++i)
    for (int j=0 ; j<m.cols(); ++j)
    {
      gsl_matrix_complex_set(res, i, j,
        gsl_complex_rect(m(i,j).real(), m(i,j).imag()));
    }
}

template<typename MatrixType>
void convert(const gsl_matrix_complex* m, MatrixType& res)
{
  res.resize(int(m->size1), int(m->size2));
  for (int i=0 ; i<res.rows() ; ++i)
    for (int j=0 ; j<res.cols(); ++j)
      res(i,j) = typename MatrixType::Scalar(
        GSL_REAL(gsl_matrix_complex_get(m,i,j)),
        GSL_IMAG(gsl_matrix_complex_get(m,i,j)));
}

template<typename VectorType>
void convert(const VectorType& m, gsl_vector_complex* &res)
{
  res = gsl_vector_complex_alloc(m.size());
  for (int i=0 ; i<m.size() ; ++i)
      gsl_vector_complex_set(res, i, gsl_complex_rect(m[i].real(), m[i].imag()));
}

template<typename VectorType>
void convert(const gsl_vector_complex* m, VectorType& res)
{
  res.resize(m->size);
  for (int i=0 ; i<res.rows() ; ++i)
    res[i] = typename VectorType::Scalar(
        GSL_REAL(gsl_vector_complex_get(m, i)),
        GSL_IMAG(gsl_vector_complex_get(m, i)));
}

}

#endif // EIGEN_GSL_HELPER
