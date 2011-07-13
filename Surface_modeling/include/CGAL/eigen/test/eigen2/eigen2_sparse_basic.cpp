// This file is part of Eigen, a lightweight C++ template library
// for linear algebra. Eigen itself is part of the KDE project.
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

template<typename SetterType,typename DenseType, typename Scalar, int Options>
bool test_random_setter(SparseMatrix<Scalar,Options>& sm, const DenseType& ref, const std::vector<Vector2i>& nonzeroCoords)
{
  typedef SparseMatrix<Scalar,Options> SparseType;
  {
    sm.setZero();
    SetterType w(sm);
    std::vector<Vector2i> remaining = nonzeroCoords;
    while(!remaining.empty())
    {
      int i = ei_random<int>(0,remaining.size()-1);
      w(remaining[i].x(),remaining[i].y()) = ref.coeff(remaining[i].x(),remaining[i].y());
      remaining[i] = remaining.back();
      remaining.pop_back();
    }
  }
  return sm.isApprox(ref);
}

template<typename SetterType,typename DenseType, typename T>
bool test_random_setter(DynamicSparseMatrix<T>& sm, const DenseType& ref, const std::vector<Vector2i>& nonzeroCoords)
{
  sm.setZero();
  std::vector<Vector2i> remaining = nonzeroCoords;
  while(!remaining.empty())
  {
    int i = ei_random<int>(0,remaining.size()-1);
    sm.coeffRef(remaining[i].x(),remaining[i].y()) = ref.coeff(remaining[i].x(),remaining[i].y());
    remaining[i] = remaining.back();
    remaining.pop_back();
  }
  return sm.isApprox(ref);
}

template<typename SparseMatrixType> void sparse_basic(const SparseMatrixType& ref)
{
  const int rows = ref.rows();
  const int cols = ref.cols();
  typedef typename SparseMatrixType::Scalar Scalar;
  enum { Flags = SparseMatrixType::Flags };
  
  double density = std::max(8./(rows*cols), 0.01);
  typedef Matrix<Scalar,Dynamic,Dynamic> DenseMatrix;
  typedef Matrix<Scalar,Dynamic,1> DenseVector;
  Scalar eps = 1e-6;

  SparseMatrixType m(rows, cols);
  DenseMatrix refMat = DenseMatrix::Zero(rows, cols);
  DenseVector vec1 = DenseVector::Random(rows);
  Scalar s1 = ei_random<Scalar>();

  std::vector<Vector2i> zeroCoords;
  std::vector<Vector2i> nonzeroCoords;
  initSparse<Scalar>(density, refMat, m, 0, &zeroCoords, &nonzeroCoords);
  
  if (zeroCoords.size()==0 || nonzeroCoords.size()==0)
    return;

  // test coeff and coeffRef
  for (int i=0; i<(int)zeroCoords.size(); ++i)
  {
    VERIFY_IS_MUCH_SMALLER_THAN( m.coeff(zeroCoords[i].x(),zeroCoords[i].y()), eps );
    if(ei_is_same_type<SparseMatrixType,SparseMatrix<Scalar,Flags> >::ret)
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
    int j = ei_random<int>(0,cols-1);
    int i = ei_random<int>(0,rows-1);
    int w = ei_random<int>(1,cols-j-1);
    int h = ei_random<int>(1,rows-i-1);

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

  // test SparseSetters
  // coherent setter
  // TODO extend the MatrixSetter
//   {
//     m.setZero();
//     VERIFY_IS_NOT_APPROX(m, refMat);
//     SparseSetter<SparseMatrixType, FullyCoherentAccessPattern> w(m);
//     for (int i=0; i<nonzeroCoords.size(); ++i)
//     {
//       w->coeffRef(nonzeroCoords[i].x(),nonzeroCoords[i].y()) = refMat.coeff(nonzeroCoords[i].x(),nonzeroCoords[i].y());
//     }
//   }
//   VERIFY_IS_APPROX(m, refMat);

  // random setter
//   {
//     m.setZero();
//     VERIFY_IS_NOT_APPROX(m, refMat);
//     SparseSetter<SparseMatrixType, RandomAccessPattern> w(m);
//     std::vector<Vector2i> remaining = nonzeroCoords;
//     while(!remaining.empty())
//     {
//       int i = ei_random<int>(0,remaining.size()-1);
//       w->coeffRef(remaining[i].x(),remaining[i].y()) = refMat.coeff(remaining[i].x(),remaining[i].y());
//       remaining[i] = remaining.back();
//       remaining.pop_back();
//     }
//   }
//   VERIFY_IS_APPROX(m, refMat);

    VERIFY(( test_random_setter<RandomSetter<SparseMatrixType, StdMapTraits> >(m,refMat,nonzeroCoords) ));
    #ifdef EIGEN_UNORDERED_MAP_SUPPORT
    VERIFY(( test_random_setter<RandomSetter<SparseMatrixType, StdUnorderedMapTraits> >(m,refMat,nonzeroCoords) ));
    #endif
    #ifdef _DENSE_HASH_MAP_H_
    VERIFY(( test_random_setter<RandomSetter<SparseMatrixType, GoogleDenseHashMapTraits> >(m,refMat,nonzeroCoords) ));
    #endif
    #ifdef _SPARSE_HASH_MAP_H_
    VERIFY(( test_random_setter<RandomSetter<SparseMatrixType, GoogleSparseHashMapTraits> >(m,refMat,nonzeroCoords) ));
    #endif

    // test fillrand
    {
      DenseMatrix m1(rows,cols);
      m1.setZero();
      SparseMatrixType m2(rows,cols);
      m2.startFill();
      for (int j=0; j<cols; ++j)
      {
        for (int k=0; k<rows/2; ++k)
        {
          int i = ei_random<int>(0,rows-1);
          if (m1.coeff(i,j)==Scalar(0))
            m2.fillrand(i,j) = m1(i,j) = ei_random<Scalar>();
        }
      }
      m2.endFill();
      VERIFY_IS_APPROX(m2,m1);
    }
  
  // test RandomSetter
  /*{
    SparseMatrixType m1(rows,cols), m2(rows,cols);
    DenseMatrix refM1 = DenseMatrix::Zero(rows, rows);
    initSparse<Scalar>(density, refM1, m1);
    {
      Eigen::RandomSetter<SparseMatrixType > setter(m2);
      for (int j=0; j<m1.outerSize(); ++j)
        for (typename SparseMatrixType::InnerIterator i(m1,j); i; ++i)
          setter(i.index(), j) = i.value();
    }
    VERIFY_IS_APPROX(m1, m2);
  }*/
//   std::cerr << m.transpose() << "\n\n"  << refMat.transpose() << "\n\n";
//   VERIFY_IS_APPROX(m, refMat);

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
    VERIFY_IS_APPROX(m3.cwise()*(m1+m2), refM3.cwise()*(refM1+refM2));
    VERIFY_IS_APPROX(m1*s1-m2, refM1*s1-refM2);

    VERIFY_IS_APPROX(m1*=s1, refM1*=s1);
    VERIFY_IS_APPROX(m1/=s1, refM1/=s1);
    
    VERIFY_IS_APPROX(m1+=m2, refM1+=refM2);
    VERIFY_IS_APPROX(m1-=m2, refM1-=refM2);
    
    VERIFY_IS_APPROX(m1.col(0).eigen2_dot(refM2.row(0)), refM1.col(0).eigen2_dot(refM2.row(0)));
    
    refM4.setRandom();
    // sparse cwise* dense
    VERIFY_IS_APPROX(m3.cwise()*refM4, refM3.cwise()*refM4);
//     VERIFY_IS_APPROX(m3.cwise()/refM4, refM3.cwise()/refM4);
  }

  // test innerVector()
  {
    DenseMatrix refMat2 = DenseMatrix::Zero(rows, rows);
    SparseMatrixType m2(rows, rows);
    initSparse<Scalar>(density, refMat2, m2);
    int j0 = ei_random(0,rows-1);
    int j1 = ei_random(0,rows-1);
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
    int j0 = ei_random(0,rows-2);
    int j1 = ei_random(0,rows-2);
    int n0 = ei_random<int>(1,rows-std::max(j0,j1));
    VERIFY_IS_APPROX(m2.innerVectors(j0,n0), refMat2.block(0,j0,rows,n0));
    VERIFY_IS_APPROX(m2.innerVectors(j0,n0)+m2.innerVectors(j1,n0),
                     refMat2.block(0,j0,rows,n0)+refMat2.block(0,j1,rows,n0));
    //m2.innerVectors(j0,n0) = m2.innerVectors(j0,n0) + m2.innerVectors(j1,n0);
    //refMat2.block(0,j0,rows,n0) = refMat2.block(0,j0,rows,n0) + refMat2.block(0,j1,rows,n0);
  }

  // test transpose
  {
    DenseMatrix refMat2 = DenseMatrix::Zero(rows, rows);
    SparseMatrixType m2(rows, rows);
    initSparse<Scalar>(density, refMat2, m2);
    VERIFY_IS_APPROX(m2.transpose().eval(), refMat2.transpose().eval());
    VERIFY_IS_APPROX(m2.transpose(), refMat2.transpose());
  }
  
  // test prune
  {
    SparseMatrixType m2(rows, rows);
    DenseMatrix refM2(rows, rows);
    refM2.setZero();
    int countFalseNonZero = 0;
    int countTrueNonZero = 0;
    m2.startFill();
    for (int j=0; j<m2.outerSize(); ++j)
      for (int i=0; i<m2.innerSize(); ++i)
      {
        float x = ei_random<float>(0,1);
        if (x<0.1)
        {
          // do nothing
        }
        else if (x<0.5)
        {
          countFalseNonZero++;
          m2.fill(i,j) = Scalar(0);
        }
        else
        {
          countTrueNonZero++;
          m2.fill(i,j) = refM2(i,j) = Scalar(1);
        }
      }
    m2.endFill();
    VERIFY(countFalseNonZero+countTrueNonZero == m2.nonZeros());
    VERIFY_IS_APPROX(m2, refM2);
    m2.prune(1);
    VERIFY(countTrueNonZero==m2.nonZeros());
    VERIFY_IS_APPROX(m2, refM2);
  }
}

void test_eigen2_sparse_basic()
{
  for(int i = 0; i < g_repeat; i++) {
    CALL_SUBTEST_1( sparse_basic(SparseMatrix<double>(8, 8)) );
    CALL_SUBTEST_2( sparse_basic(SparseMatrix<std::complex<double> >(16, 16)) );
    CALL_SUBTEST_1( sparse_basic(SparseMatrix<double>(33, 33)) );
    
    CALL_SUBTEST_3( sparse_basic(DynamicSparseMatrix<double>(8, 8)) );
  }
}
