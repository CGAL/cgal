// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2008 Gael Guennebaud <gael.guennebaud@inria.fr>
// Copyright (C) 2009 Benoit Jacob <jacob.benoit.1@gmail.com>
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

// discard stack allocation as that too bypasses malloc
#define EIGEN_STACK_ALLOCATION_LIMIT 0
#define EIGEN_RUNTIME_NO_MALLOC
#include "main.h"
#include <Eigen/SVD>

template<typename MatrixType, int QRPreconditioner>
void jacobisvd_check_full(const MatrixType& m, const JacobiSVD<MatrixType, QRPreconditioner>& svd)
{
  typedef typename MatrixType::Index Index;
  Index rows = m.rows();
  Index cols = m.cols();

  enum {
    RowsAtCompileTime = MatrixType::RowsAtCompileTime,
    ColsAtCompileTime = MatrixType::ColsAtCompileTime
  };

  typedef typename MatrixType::Scalar Scalar;
  typedef typename NumTraits<Scalar>::Real RealScalar;
  typedef Matrix<Scalar, RowsAtCompileTime, RowsAtCompileTime> MatrixUType;
  typedef Matrix<Scalar, ColsAtCompileTime, ColsAtCompileTime> MatrixVType;
  typedef Matrix<Scalar, RowsAtCompileTime, 1> ColVectorType;
  typedef Matrix<Scalar, ColsAtCompileTime, 1> InputVectorType;

  MatrixType sigma = MatrixType::Zero(rows,cols);
  sigma.diagonal() = svd.singularValues().template cast<Scalar>();
  MatrixUType u = svd.matrixU();
  MatrixVType v = svd.matrixV();

  VERIFY_IS_APPROX(m, u * sigma * v.adjoint());
  VERIFY_IS_UNITARY(u);
  VERIFY_IS_UNITARY(v);
}

template<typename MatrixType, int QRPreconditioner>
void jacobisvd_compare_to_full(const MatrixType& m,
                               unsigned int computationOptions,
                               const JacobiSVD<MatrixType, QRPreconditioner>& referenceSvd)
{
  typedef typename MatrixType::Index Index;
  Index rows = m.rows();
  Index cols = m.cols();
  Index diagSize = std::min(rows, cols);

  JacobiSVD<MatrixType, QRPreconditioner> svd(m, computationOptions);

  VERIFY_IS_APPROX(svd.singularValues(), referenceSvd.singularValues());
  if(computationOptions & ComputeFullU)
    VERIFY_IS_APPROX(svd.matrixU(), referenceSvd.matrixU());
  if(computationOptions & ComputeThinU)
    VERIFY_IS_APPROX(svd.matrixU(), referenceSvd.matrixU().leftCols(diagSize));
  if(computationOptions & ComputeFullV)
    VERIFY_IS_APPROX(svd.matrixV(), referenceSvd.matrixV());
  if(computationOptions & ComputeThinV)
    VERIFY_IS_APPROX(svd.matrixV(), referenceSvd.matrixV().leftCols(diagSize));
}

template<typename MatrixType, int QRPreconditioner>
void jacobisvd_solve(const MatrixType& m, unsigned int computationOptions)
{
  typedef typename MatrixType::Scalar Scalar;
  typedef typename MatrixType::Index Index;
  Index rows = m.rows();
  Index cols = m.cols();

  enum {
    RowsAtCompileTime = MatrixType::RowsAtCompileTime,
    ColsAtCompileTime = MatrixType::ColsAtCompileTime
  };

  typedef Matrix<Scalar, RowsAtCompileTime, Dynamic> RhsType;
  typedef Matrix<Scalar, ColsAtCompileTime, Dynamic> SolutionType;

  RhsType rhs = RhsType::Random(rows, internal::random<Index>(1, cols));
  JacobiSVD<MatrixType, QRPreconditioner> svd(m, computationOptions);
  SolutionType x = svd.solve(rhs);
  // evaluate normal equation which works also for least-squares solutions
  VERIFY_IS_APPROX(m.adjoint()*m*x,m.adjoint()*rhs);
}

template<typename MatrixType, int QRPreconditioner>
void jacobisvd_test_all_computation_options(const MatrixType& m)
{
  if (QRPreconditioner == NoQRPreconditioner && m.rows() != m.cols())
    return;
  JacobiSVD<MatrixType, QRPreconditioner> fullSvd(m, ComputeFullU|ComputeFullV);

  jacobisvd_check_full(m, fullSvd);
  jacobisvd_solve<MatrixType, QRPreconditioner>(m, ComputeFullU | ComputeFullV);

  if(QRPreconditioner == FullPivHouseholderQRPreconditioner)
    return;

  jacobisvd_compare_to_full(m, ComputeFullU, fullSvd);
  jacobisvd_compare_to_full(m, ComputeFullV, fullSvd);
  jacobisvd_compare_to_full(m, 0, fullSvd);

  if (MatrixType::ColsAtCompileTime == Dynamic) {
    // thin U/V are only available with dynamic number of columns
    jacobisvd_compare_to_full(m, ComputeFullU|ComputeThinV, fullSvd);
    jacobisvd_compare_to_full(m,              ComputeThinV, fullSvd);
    jacobisvd_compare_to_full(m, ComputeThinU|ComputeFullV, fullSvd);
    jacobisvd_compare_to_full(m, ComputeThinU             , fullSvd);
    jacobisvd_compare_to_full(m, ComputeThinU|ComputeThinV, fullSvd);
    jacobisvd_solve<MatrixType, QRPreconditioner>(m, ComputeFullU | ComputeThinV);
    jacobisvd_solve<MatrixType, QRPreconditioner>(m, ComputeThinU | ComputeFullV);
    jacobisvd_solve<MatrixType, QRPreconditioner>(m, ComputeThinU | ComputeThinV);
  }
}

template<typename MatrixType>
void jacobisvd(const MatrixType& a = MatrixType(), bool pickrandom = true)
{
  MatrixType m = pickrandom ? MatrixType::Random(a.rows(), a.cols()) : a;

  jacobisvd_test_all_computation_options<MatrixType, FullPivHouseholderQRPreconditioner>(m);
  jacobisvd_test_all_computation_options<MatrixType, ColPivHouseholderQRPreconditioner>(m);
  jacobisvd_test_all_computation_options<MatrixType, HouseholderQRPreconditioner>(m);
  jacobisvd_test_all_computation_options<MatrixType, NoQRPreconditioner>(m);
}

template<typename MatrixType> void jacobisvd_verify_assert(const MatrixType& m)
{
  typedef typename MatrixType::Scalar Scalar;
  typedef typename MatrixType::Index Index;
  Index rows = m.rows();
  Index cols = m.cols();

  enum {
    RowsAtCompileTime = MatrixType::RowsAtCompileTime,
    ColsAtCompileTime = MatrixType::ColsAtCompileTime
  };

  typedef Matrix<Scalar, RowsAtCompileTime, 1> RhsType;

  RhsType rhs(rows);

  JacobiSVD<MatrixType> svd;
  VERIFY_RAISES_ASSERT(svd.matrixU())
  VERIFY_RAISES_ASSERT(svd.singularValues())
  VERIFY_RAISES_ASSERT(svd.matrixV())
  VERIFY_RAISES_ASSERT(svd.solve(rhs))

  MatrixType a = MatrixType::Zero(rows, cols);
  a.setZero();
  svd.compute(a, 0);
  VERIFY_RAISES_ASSERT(svd.matrixU())
  VERIFY_RAISES_ASSERT(svd.matrixV())
  svd.singularValues();
  VERIFY_RAISES_ASSERT(svd.solve(rhs))

  if (ColsAtCompileTime == Dynamic)
  {
    svd.compute(a, ComputeThinU);
    svd.matrixU();
    VERIFY_RAISES_ASSERT(svd.matrixV())
    VERIFY_RAISES_ASSERT(svd.solve(rhs))

    svd.compute(a, ComputeThinV);
    svd.matrixV();
    VERIFY_RAISES_ASSERT(svd.matrixU())
    VERIFY_RAISES_ASSERT(svd.solve(rhs))

    JacobiSVD<MatrixType, FullPivHouseholderQRPreconditioner> svd_fullqr;
    VERIFY_RAISES_ASSERT(svd_fullqr.compute(a, ComputeFullU|ComputeThinV))
    VERIFY_RAISES_ASSERT(svd_fullqr.compute(a, ComputeThinU|ComputeThinV))
    VERIFY_RAISES_ASSERT(svd_fullqr.compute(a, ComputeThinU|ComputeFullV))
  }
  else
  {
    VERIFY_RAISES_ASSERT(svd.compute(a, ComputeThinU))
    VERIFY_RAISES_ASSERT(svd.compute(a, ComputeThinV))
  }
}

template<typename MatrixType>
void jacobisvd_method()
{
  enum { Size = MatrixType::RowsAtCompileTime };
  typedef typename MatrixType::RealScalar RealScalar;
  typedef Matrix<RealScalar, Size, 1> RealVecType;
  MatrixType m = MatrixType::Identity();
  VERIFY_IS_APPROX(m.jacobiSvd().singularValues(), RealVecType::Ones());
  VERIFY_RAISES_ASSERT(m.jacobiSvd().matrixU());
  VERIFY_RAISES_ASSERT(m.jacobiSvd().matrixV());
  VERIFY_IS_APPROX(m.jacobiSvd(ComputeFullU|ComputeFullV).solve(m), m);
}

// work around stupid msvc error when constructing at compile time an expression that involves
// a division by zero, even if the numeric type has floating point
template<typename Scalar>
EIGEN_DONT_INLINE Scalar zero() { return Scalar(0); }

// workaround aggressive optimization in ICC
template<typename T> EIGEN_DONT_INLINE  T sub(T a, T b) { return a - b; }

template<typename MatrixType>
void jacobisvd_inf_nan()
{
  // all this function does is verify we don't iterate infinitely on nan/inf values

  JacobiSVD<MatrixType> svd;
  typedef typename MatrixType::Scalar Scalar;
  Scalar some_inf = Scalar(1) / zero<Scalar>();
  VERIFY(sub(some_inf, some_inf) != sub(some_inf, some_inf));
  svd.compute(MatrixType::Constant(10,10,some_inf), ComputeFullU | ComputeFullV);

  Scalar some_nan = zero<Scalar>() / zero<Scalar>();
  VERIFY(some_nan != some_nan);
  svd.compute(MatrixType::Constant(10,10,some_nan), ComputeFullU | ComputeFullV);

  MatrixType m = MatrixType::Zero(10,10);
  m(internal::random<int>(0,9), internal::random<int>(0,9)) = some_inf;
  svd.compute(m, ComputeFullU | ComputeFullV);

  m = MatrixType::Zero(10,10);
  m(internal::random<int>(0,9), internal::random<int>(0,9)) = some_nan;
  svd.compute(m, ComputeFullU | ComputeFullV);
}

void jacobisvd_preallocate()
{
  Vector3f v(3.f, 2.f, 1.f);
  MatrixXf m = v.asDiagonal();

  internal::set_is_malloc_allowed(false);
  VERIFY_RAISES_ASSERT(VectorXf v(10);)
  JacobiSVD<MatrixXf> svd;
  internal::set_is_malloc_allowed(true);
  svd.compute(m);
  VERIFY_IS_APPROX(svd.singularValues(), v);

  JacobiSVD<MatrixXf> svd2(3,3);
  internal::set_is_malloc_allowed(false);
  svd2.compute(m);
  internal::set_is_malloc_allowed(true);
  VERIFY_IS_APPROX(svd2.singularValues(), v);
  VERIFY_RAISES_ASSERT(svd2.matrixU());
  VERIFY_RAISES_ASSERT(svd2.matrixV());
  svd2.compute(m, ComputeFullU | ComputeFullV);
  VERIFY_IS_APPROX(svd2.matrixU(), Matrix3f::Identity());
  VERIFY_IS_APPROX(svd2.matrixV(), Matrix3f::Identity());
  internal::set_is_malloc_allowed(false);
  svd2.compute(m);
  internal::set_is_malloc_allowed(true);

  JacobiSVD<MatrixXf> svd3(3,3,ComputeFullU|ComputeFullV);
  internal::set_is_malloc_allowed(false);
  svd2.compute(m);
  internal::set_is_malloc_allowed(true);
  VERIFY_IS_APPROX(svd2.singularValues(), v);
  VERIFY_IS_APPROX(svd2.matrixU(), Matrix3f::Identity());
  VERIFY_IS_APPROX(svd2.matrixV(), Matrix3f::Identity());
  internal::set_is_malloc_allowed(false);
  svd2.compute(m, ComputeFullU|ComputeFullV);
  internal::set_is_malloc_allowed(true);

  
}

void test_jacobisvd()
{
  CALL_SUBTEST_3(( jacobisvd_verify_assert(Matrix3f()) ));
  CALL_SUBTEST_4(( jacobisvd_verify_assert(Matrix4d()) ));
  CALL_SUBTEST_7(( jacobisvd_verify_assert(MatrixXf(10,12)) ));
  CALL_SUBTEST_8(( jacobisvd_verify_assert(MatrixXcd(7,5)) ));

  for(int i = 0; i < g_repeat; i++) {
    Matrix2cd m;
    m << 0, 1,
         0, 1;
    CALL_SUBTEST_1(( jacobisvd(m, false) ));
    m << 1, 0,
         1, 0;
    CALL_SUBTEST_1(( jacobisvd(m, false) ));

    Matrix2d n;
    n << 0, 0,
         0, 0;
    CALL_SUBTEST_2(( jacobisvd(n, false) ));
    n << 0, 0,
         0, 1;
    CALL_SUBTEST_2(( jacobisvd(n, false) ));
    
    CALL_SUBTEST_3(( jacobisvd<Matrix3f>() ));
    CALL_SUBTEST_4(( jacobisvd<Matrix4d>() ));
    CALL_SUBTEST_5(( jacobisvd<Matrix<float,3,5> >() ));
    CALL_SUBTEST_6(( jacobisvd<Matrix<double,Dynamic,2> >(Matrix<double,Dynamic,2>(10,2)) ));

    int r = internal::random<int>(1, 30),
        c = internal::random<int>(1, 30);
    CALL_SUBTEST_7(( jacobisvd<MatrixXf>(MatrixXf(r,c)) ));
    CALL_SUBTEST_8(( jacobisvd<MatrixXcd>(MatrixXcd(r,c)) ));
    (void) r;
    (void) c;

    // Test on inf/nan matrix
    CALL_SUBTEST_7( jacobisvd_inf_nan<MatrixXf>() );
  }

  CALL_SUBTEST_7(( jacobisvd<MatrixXf>(MatrixXf(internal::random<int>(100, 150), internal::random<int>(100, 150))) ));
  CALL_SUBTEST_8(( jacobisvd<MatrixXcd>(MatrixXcd(internal::random<int>(80, 100), internal::random<int>(80, 100))) ));

  // test matrixbase method
  CALL_SUBTEST_1(( jacobisvd_method<Matrix2cd>() ));
  CALL_SUBTEST_3(( jacobisvd_method<Matrix3f>() ));

  // Test problem size constructors
  CALL_SUBTEST_7( JacobiSVD<MatrixXf>(10,10) );

  // Check that preallocation avoids subsequent mallocs
  CALL_SUBTEST_9( jacobisvd_preallocate() );
}
