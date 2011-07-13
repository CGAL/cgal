// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2008-2009 Benoit Jacob <jacob.benoit.1@gmail.com>
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

#include "main.h"
#include <Eigen/LU>
using namespace std;

template<typename MatrixType> void lu_non_invertible()
{
  typedef typename MatrixType::Index Index;
  typedef typename MatrixType::Scalar Scalar;
  typedef typename MatrixType::RealScalar RealScalar;
  /* this test covers the following files:
     LU.h
  */
  Index rows, cols, cols2;
  if(MatrixType::RowsAtCompileTime==Dynamic)
  {
    rows = internal::random<Index>(2,200);
  }
  else
  {
    rows = MatrixType::RowsAtCompileTime;
  }
  if(MatrixType::ColsAtCompileTime==Dynamic)
  {
    cols = internal::random<Index>(2,200);
    cols2 = internal::random<int>(2,200);
  }
  else
  {
    cols2 = cols = MatrixType::ColsAtCompileTime;
  }

  enum {
    RowsAtCompileTime = MatrixType::RowsAtCompileTime,
    ColsAtCompileTime = MatrixType::ColsAtCompileTime
  };
  typedef typename internal::kernel_retval_base<FullPivLU<MatrixType> >::ReturnType KernelMatrixType;
  typedef typename internal::image_retval_base<FullPivLU<MatrixType> >::ReturnType ImageMatrixType;
  typedef Matrix<typename MatrixType::Scalar, ColsAtCompileTime, ColsAtCompileTime>
          CMatrixType;
  typedef Matrix<typename MatrixType::Scalar, RowsAtCompileTime, RowsAtCompileTime>
          RMatrixType;

  Index rank = internal::random<Index>(1, std::min(rows, cols)-1);

  // The image of the zero matrix should consist of a single (zero) column vector
  VERIFY((MatrixType::Zero(rows,cols).fullPivLu().image(MatrixType::Zero(rows,cols)).cols() == 1));

  MatrixType m1(rows, cols), m3(rows, cols2);
  CMatrixType m2(cols, cols2);
  createRandomPIMatrixOfRank(rank, rows, cols, m1);

  FullPivLU<MatrixType> lu;

  // The special value 0.01 below works well in tests. Keep in mind that we're only computing the rank
  // of singular values are either 0 or 1.
  // So it's not clear at all that the epsilon should play any role there.
  lu.setThreshold(RealScalar(0.01));
  lu.compute(m1);

  MatrixType u(rows,cols);
  u = lu.matrixLU().template triangularView<Upper>();
  RMatrixType l = RMatrixType::Identity(rows,rows);
  l.block(0,0,rows,std::min(rows,cols)).template triangularView<StrictlyLower>()
    = lu.matrixLU().block(0,0,rows,std::min(rows,cols));

  VERIFY_IS_APPROX(lu.permutationP() * m1 * lu.permutationQ(), l*u);

  KernelMatrixType m1kernel = lu.kernel();
  ImageMatrixType m1image = lu.image(m1);

  VERIFY_IS_APPROX(m1, lu.reconstructedMatrix());
  VERIFY(rank == lu.rank());
  VERIFY(cols - lu.rank() == lu.dimensionOfKernel());
  VERIFY(!lu.isInjective());
  VERIFY(!lu.isInvertible());
  VERIFY(!lu.isSurjective());
  VERIFY((m1 * m1kernel).isMuchSmallerThan(m1));
  VERIFY(m1image.fullPivLu().rank() == rank);
  VERIFY_IS_APPROX(m1 * m1.adjoint() * m1image, m1image);

  m2 = CMatrixType::Random(cols,cols2);
  m3 = m1*m2;
  m2 = CMatrixType::Random(cols,cols2);
  // test that the code, which does resize(), may be applied to an xpr
  m2.block(0,0,m2.rows(),m2.cols()) = lu.solve(m3);
  VERIFY_IS_APPROX(m3, m1*m2);
}

template<typename MatrixType> void lu_invertible()
{
  /* this test covers the following files:
     LU.h
  */
  typedef typename MatrixType::Scalar Scalar;
  typedef typename NumTraits<typename MatrixType::Scalar>::Real RealScalar;
  int size = internal::random<int>(1,200);

  MatrixType m1(size, size), m2(size, size), m3(size, size);
  FullPivLU<MatrixType> lu;
  lu.setThreshold(RealScalar(0.01));
  do {
    m1 = MatrixType::Random(size,size);
    lu.compute(m1);
  } while(!lu.isInvertible());

  VERIFY_IS_APPROX(m1, lu.reconstructedMatrix());
  VERIFY(0 == lu.dimensionOfKernel());
  VERIFY(lu.kernel().cols() == 1); // the kernel() should consist of a single (zero) column vector
  VERIFY(size == lu.rank());
  VERIFY(lu.isInjective());
  VERIFY(lu.isSurjective());
  VERIFY(lu.isInvertible());
  VERIFY(lu.image(m1).fullPivLu().isInvertible());
  m3 = MatrixType::Random(size,size);
  m2 = lu.solve(m3);
  VERIFY_IS_APPROX(m3, m1*m2);
  VERIFY_IS_APPROX(m2, lu.inverse()*m3);
}

template<typename MatrixType> void lu_partial_piv()
{
  /* this test covers the following files:
     PartialPivLU.h
  */
  typedef typename MatrixType::Index Index;
  typedef typename MatrixType::Scalar Scalar;
  typedef typename NumTraits<typename MatrixType::Scalar>::Real RealScalar;
  Index rows = internal::random<Index>(1,4);
  Index cols = rows;

  MatrixType m1(cols, rows);
  m1.setRandom();
  PartialPivLU<MatrixType> plu(m1);

  VERIFY_IS_APPROX(m1, plu.reconstructedMatrix());
}

template<typename MatrixType> void lu_verify_assert()
{
  MatrixType tmp;

  FullPivLU<MatrixType> lu;
  VERIFY_RAISES_ASSERT(lu.matrixLU())
  VERIFY_RAISES_ASSERT(lu.permutationP())
  VERIFY_RAISES_ASSERT(lu.permutationQ())
  VERIFY_RAISES_ASSERT(lu.kernel())
  VERIFY_RAISES_ASSERT(lu.image(tmp))
  VERIFY_RAISES_ASSERT(lu.solve(tmp))
  VERIFY_RAISES_ASSERT(lu.determinant())
  VERIFY_RAISES_ASSERT(lu.rank())
  VERIFY_RAISES_ASSERT(lu.dimensionOfKernel())
  VERIFY_RAISES_ASSERT(lu.isInjective())
  VERIFY_RAISES_ASSERT(lu.isSurjective())
  VERIFY_RAISES_ASSERT(lu.isInvertible())
  VERIFY_RAISES_ASSERT(lu.inverse())

  PartialPivLU<MatrixType> plu;
  VERIFY_RAISES_ASSERT(plu.matrixLU())
  VERIFY_RAISES_ASSERT(plu.permutationP())
  VERIFY_RAISES_ASSERT(plu.solve(tmp))
  VERIFY_RAISES_ASSERT(plu.determinant())
  VERIFY_RAISES_ASSERT(plu.inverse())
}

void test_lu()
{
  for(int i = 0; i < g_repeat; i++) {
    CALL_SUBTEST_1( lu_non_invertible<Matrix3f>() );
    CALL_SUBTEST_1( lu_verify_assert<Matrix3f>() );

    CALL_SUBTEST_2( (lu_non_invertible<Matrix<double, 4, 6> >()) );
    CALL_SUBTEST_2( (lu_verify_assert<Matrix<double, 4, 6> >()) );

    CALL_SUBTEST_3( lu_non_invertible<MatrixXf>() );
    CALL_SUBTEST_3( lu_invertible<MatrixXf>() );
    CALL_SUBTEST_3( lu_verify_assert<MatrixXf>() );

    CALL_SUBTEST_4( lu_non_invertible<MatrixXd>() );
    CALL_SUBTEST_4( lu_invertible<MatrixXd>() );
    CALL_SUBTEST_4( lu_partial_piv<MatrixXd>() );
    CALL_SUBTEST_4( lu_verify_assert<MatrixXd>() );

    CALL_SUBTEST_5( lu_non_invertible<MatrixXcf>() );
    CALL_SUBTEST_5( lu_invertible<MatrixXcf>() );
    CALL_SUBTEST_5( lu_verify_assert<MatrixXcf>() );

    CALL_SUBTEST_6( lu_non_invertible<MatrixXcd>() );
    CALL_SUBTEST_6( lu_invertible<MatrixXcd>() );
    CALL_SUBTEST_6( lu_partial_piv<MatrixXcd>() );
    CALL_SUBTEST_6( lu_verify_assert<MatrixXcd>() );

    CALL_SUBTEST_7(( lu_non_invertible<Matrix<float,Dynamic,16> >() ));

    // Test problem size constructors
    CALL_SUBTEST_9( PartialPivLU<MatrixXf>(10) );
    CALL_SUBTEST_9( FullPivLU<MatrixXf>(10, 20); );
  }
}
