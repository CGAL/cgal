// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2008 Daniel Gomez Ferro <dgomezferro@gmail.com>
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

#include "sparse.h"

template<typename SparseMatrixType> void sparse_basic(const SparseMatrixType& ref)
{
  typedef typename SparseMatrixType::Index Index;

  const Index rows = ref.rows();
  const Index cols = ref.cols();
  typedef typename SparseMatrixType::Scalar Scalar;
  enum { Flags = SparseMatrixType::Flags };

  double density = std::max(8./(rows*cols), 0.01);
  typedef Matrix<Scalar,Dynamic,Dynamic> DenseMatrix;
  typedef Matrix<Scalar,Dynamic,1> DenseVector;
  Scalar eps = 1e-6;

  SparseMatrixType m(rows, cols);
  DenseMatrix refMat = DenseMatrix::Zero(rows, cols);
  DenseVector vec1 = DenseVector::Random(rows);
  Scalar s1 = internal::random<Scalar>();

  std::vector<Vector2i> zeroCoords;
  std::vector<Vector2i> nonzeroCoords;
  initSparse<Scalar>(density, refMat, m, 0, &zeroCoords, &nonzeroCoords);

  if (zeroCoords.size()==0 || nonzeroCoords.size()==0)
    return;

  // test coeff and coeffRef
  for (int i=0; i<(int)zeroCoords.size(); ++i)
  {
    VERIFY_IS_MUCH_SMALLER_THAN( m.coeff(zeroCoords[i].x(),zeroCoords[i].y()), eps );
    if(internal::is_same<SparseMatrixType,SparseMatrix<Scalar,Flags> >::value)
      VERIFY_RAISES_ASSERT( m.coeffRef(zeroCoords[0].x(),zeroCoords[0].y()) = 5 );
  }
  VERIFY_IS_APPROX(m, refMat);

  m.coeffRef(nonzeroCoords[0].x(), nonzeroCoords[0].y()) = Scalar(5);
  refMat.coeffRef(nonzeroCoords[0].x(), nonzeroCoords[0].y()) = Scalar(5);

  VERIFY_IS_APPROX(m, refMat);
  /*
  // test InnerIterators and Block expressions
  for (int t=0; t<10; ++t)
  {
    int j = internal::random<int>(0,cols-1);
    int i = internal::random<int>(0,rows-1);
    int w = internal::random<int>(1,cols-j-1);
    int h = internal::random<int>(1,rows-i-1);

//     VERIFY_IS_APPROX(m.block(i,j,h,w), refMat.block(i,j,h,w));
    for(int c=0; c<w; c++)
    {
      VERIFY_IS_APPROX(m.block(i,j,h,w).col(c), refMat.block(i,j,h,w).col(c));
      for(int r=0; r<h; r++)
      {
//         VERIFY_IS_APPROX(m.block(i,j,h,w).col(c).coeff(r), refMat.block(i,j,h,w).col(c).coeff(r));
      }
    }
//     for(int r=0; r<h; r++)
//     {
//       VERIFY_IS_APPROX(m.block(i,j,h,w).row(r), refMat.block(i,j,h,w).row(r));
//       for(int c=0; c<w; c++)
//       {
//         VERIFY_IS_APPROX(m.block(i,j,h,w).row(r).coeff(c), refMat.block(i,j,h,w).row(r).coeff(c));
//       }
//     }
  }

  for(int c=0; c<cols; c++)
  {
    VERIFY_IS_APPROX(m.col(c) + m.col(c), (m + m).col(c));
    VERIFY_IS_APPROX(m.col(c) + m.col(c), refMat.col(c) + refMat.col(c));
  }

  for(int r=0; r<rows; r++)
  {
    VERIFY_IS_APPROX(m.row(r) + m.row(r), (m + m).row(r));
    VERIFY_IS_APPROX(m.row(r) + m.row(r), refMat.row(r) + refMat.row(r));
  }
  */

    // test insert (inner random)
    {
      DenseMatrix m1(rows,cols);
      m1.setZero();
      SparseMatrixType m2(rows,cols);
      m2.reserve(10);
      for (int j=0; j<cols; ++j)
      {
        for (int k=0; k<rows/2; ++k)
        {
          int i = internal::random<int>(0,rows-1);
          if (m1.coeff(i,j)==Scalar(0))
            m2.insert(i,j) = m1(i,j) = internal::random<Scalar>();
        }
      }
      m2.finalize();
      VERIFY_IS_APPROX(m2,m1);
    }

    // test insert (fully random)
    {
      DenseMatrix m1(rows,cols);
      m1.setZero();
      SparseMatrixType m2(rows,cols);
      m2.reserve(10);
      for (int k=0; k<rows*cols; ++k)
      {
        int i = internal::random<int>(0,rows-1);
        int j = internal::random<int>(0,cols-1);
        if (m1.coeff(i,j)==Scalar(0))
          m2.insert(i,j) = m1(i,j) = internal::random<Scalar>();
      }
      m2.finalize();
      VERIFY_IS_APPROX(m2,m1);
    }

  // test basic computations
  {
    DenseMatrix refM1 = DenseMatrix::Zero(rows, rows);
    DenseMatrix refM2 = DenseMatrix::Zero(rows, rows);
    DenseMatrix refM3 = DenseMatrix::Zero(rows, rows);
    DenseMatrix refM4 = DenseMatrix::Zero(rows, rows);
    SparseMatrixType m1(rows, rows);
    SparseMatrixType m2(rows, rows);
    SparseMatrixType m3(rows, rows);
    SparseMatrixType m4(rows, rows);
    initSparse<Scalar>(density, refM1, m1);
    initSparse<Scalar>(density, refM2, m2);
    initSparse<Scalar>(density, refM3, m3);
    initSparse<Scalar>(density, refM4, m4);

    VERIFY_IS_APPROX(m1+m2, refM1+refM2);
    VERIFY_IS_APPROX(m1+m2+m3, refM1+refM2+refM3);
    VERIFY_IS_APPROX(m3.cwiseProduct(m1+m2), refM3.cwiseProduct(refM1+refM2));
    VERIFY_IS_APPROX(m1*s1-m2, refM1*s1-refM2);

    VERIFY_IS_APPROX(m1*=s1, refM1*=s1);
    VERIFY_IS_APPROX(m1/=s1, refM1/=s1);

    VERIFY_IS_APPROX(m1+=m2, refM1+=refM2);
    VERIFY_IS_APPROX(m1-=m2, refM1-=refM2);

    VERIFY_IS_APPROX(m1.col(0).dot(refM2.row(0)), refM1.col(0).dot(refM2.row(0)));

    refM4.setRandom();
    // sparse cwise* dense
    VERIFY_IS_APPROX(m3.cwiseProduct(refM4), refM3.cwiseProduct(refM4));
//     VERIFY_IS_APPROX(m3.cwise()/refM4, refM3.cwise()/refM4);
  }

  // test transpose
  {
    DenseMatrix refMat2 = DenseMatrix::Zero(rows, rows);
    SparseMatrixType m2(rows, rows);
    initSparse<Scalar>(density, refMat2, m2);
    VERIFY_IS_APPROX(m2.transpose().eval(), refMat2.transpose().eval());
    VERIFY_IS_APPROX(m2.transpose(), refMat2.transpose());

    VERIFY_IS_APPROX(SparseMatrixType(m2.adjoint()), refMat2.adjoint());
  }

  // test innerVector()
  {
    DenseMatrix refMat2 = DenseMatrix::Zero(rows, rows);
    SparseMatrixType m2(rows, rows);
    initSparse<Scalar>(density, refMat2, m2);
    int j0 = internal::random(0,rows-1);
    int j1 = internal::random(0,rows-1);
    VERIFY_IS_APPROX(m2.innerVector(j0), refMat2.col(j0));
    VERIFY_IS_APPROX(m2.innerVector(j0)+m2.innerVector(j1), refMat2.col(j0)+refMat2.col(j1));
    //m2.innerVector(j0) = 2*m2.innerVector(j1);
    //refMat2.col(j0) = 2*refMat2.col(j1);
    //VERIFY_IS_APPROX(m2, refMat2);
  }

  // test innerVectors()
  {
    DenseMatrix refMat2 = DenseMatrix::Zero(rows, rows);
    SparseMatrixType m2(rows, rows);
    initSparse<Scalar>(density, refMat2, m2);
    int j0 = internal::random(0,rows-2);
    int j1 = internal::random(0,rows-2);
    int n0 = internal::random<int>(1,rows-std::max(j0,j1));
    VERIFY_IS_APPROX(m2.innerVectors(j0,n0), refMat2.block(0,j0,rows,n0));
    VERIFY_IS_APPROX(m2.innerVectors(j0,n0)+m2.innerVectors(j1,n0),
                     refMat2.block(0,j0,rows,n0)+refMat2.block(0,j1,rows,n0));
    //m2.innerVectors(j0,n0) = m2.innerVectors(j0,n0) + m2.innerVectors(j1,n0);
    //refMat2.block(0,j0,rows,n0) = refMat2.block(0,j0,rows,n0) + refMat2.block(0,j1,rows,n0);
  }

  // test prune
  {
    SparseMatrixType m2(rows, rows);
    DenseMatrix refM2(rows, rows);
    refM2.setZero();
    int countFalseNonZero = 0;
    int countTrueNonZero = 0;
    for (int j=0; j<m2.outerSize(); ++j)
    {
      m2.startVec(j);
      for (int i=0; i<m2.innerSize(); ++i)
      {
        float x = internal::random<float>(0,1);
        if (x<0.1)
        {
          // do nothing
        }
        else if (x<0.5)
        {
          countFalseNonZero++;
          m2.insertBackByOuterInner(j,i) = Scalar(0);
        }
        else
        {
          countTrueNonZero++;
          m2.insertBackByOuterInner(j,i) = refM2(i,j) = Scalar(1);
        }
      }
    }
    m2.finalize();
    VERIFY(countFalseNonZero+countTrueNonZero == m2.nonZeros());
    VERIFY_IS_APPROX(m2, refM2);
    m2.prune(Scalar(1));
    VERIFY(countTrueNonZero==m2.nonZeros());
    VERIFY_IS_APPROX(m2, refM2);
  }
  
  // test selfadjointView
  {
    DenseMatrix refMat2(rows, rows), refMat3(rows, rows);
    SparseMatrixType m2(rows, rows), m3(rows, rows);
    initSparse<Scalar>(density, refMat2, m2);
    refMat3 = refMat2.template selfadjointView<Lower>();
    m3 = m2.template selfadjointView<Lower>();
    VERIFY_IS_APPROX(m3, refMat3);
  }
  
  // test sparseView
  {
    DenseMatrix refMat2 = DenseMatrix::Zero(rows, rows);
    SparseMatrixType m2(rows, rows);
    initSparse<Scalar>(density, refMat2, m2);
    VERIFY_IS_APPROX(m2.eval(), refMat2.sparseView().eval());
  }
}

void test_sparse_basic()
{
  for(int i = 0; i < g_repeat; i++) {
    CALL_SUBTEST_1( sparse_basic(SparseMatrix<double>(8, 8)) );
    CALL_SUBTEST_2( sparse_basic(SparseMatrix<std::complex<double> >(16, 16)) );
    CALL_SUBTEST_1( sparse_basic(SparseMatrix<double>(33, 33)) );

    CALL_SUBTEST_3( sparse_basic(DynamicSparseMatrix<double>(8, 8)) );
  }
}
