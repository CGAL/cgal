// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2008 Gael Guennebaud <gael.guennebaud@inria.fr>
// Copyright (C) 2006-2008 Benoit Jacob <jacob.benoit.1@gmail.com>
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

// this hack is needed to make this file compiles with -pedantic (gcc)
#ifdef __GNUC__
#define throw(X)
#endif
// discard stack allocation as that too bypasses malloc
#define EIGEN_STACK_ALLOCATION_LIMIT 0
// any heap allocation will raise an assert
#define EIGEN_NO_MALLOC

#include "main.h"
#include <Eigen/Cholesky>
#include <Eigen/Eigenvalues>
#include <Eigen/LU>
#include <Eigen/QR>
#include <Eigen/SVD>

template<typename MatrixType> void nomalloc(const MatrixType& m)
{
  /* this test check no dynamic memory allocation are issued with fixed-size matrices
  */
  typedef typename MatrixType::Index Index;
  typedef typename MatrixType::Scalar Scalar;
  typedef Matrix<Scalar, MatrixType::RowsAtCompileTime, 1> VectorType;

  Index rows = m.rows();
  Index cols = m.cols();

  MatrixType m1 = MatrixType::Random(rows, cols),
             m2 = MatrixType::Random(rows, cols),
             m3(rows, cols),
             mzero = MatrixType::Zero(rows, cols),
             identity = Matrix<Scalar, MatrixType::RowsAtCompileTime, MatrixType::RowsAtCompileTime>
                              ::Identity(rows, rows),
             square = Matrix<Scalar, MatrixType::RowsAtCompileTime, MatrixType::RowsAtCompileTime>
                              ::Random(rows, rows);
  VectorType v1 = VectorType::Random(rows),
             v2 = VectorType::Random(rows),
             vzero = VectorType::Zero(rows);

  Scalar s1 = internal::random<Scalar>();

  Index r = internal::random<Index>(0, rows-1),
        c = internal::random<Index>(0, cols-1);

  VERIFY_IS_APPROX((m1+m2)*s1,              s1*m1+s1*m2);
  VERIFY_IS_APPROX((m1+m2)(r,c), (m1(r,c))+(m2(r,c)));
  VERIFY_IS_APPROX(m1.cwiseProduct(m1.block(0,0,rows,cols)), (m1.array()*m1.array()).matrix());
  VERIFY_IS_APPROX((m1*m1.transpose())*m2,  m1*(m1.transpose()*m2));
  
  m2.col(0).noalias() = m1 * m1.col(0);
  m2.col(0).noalias() -= m1.adjoint() * m1.col(0);
  m2.col(0).noalias() -= m1 * m1.row(0).adjoint();
  m2.col(0).noalias() -= m1.adjoint() * m1.row(0).adjoint();

  m2.row(0).noalias() = m1.row(0) * m1;
  m2.row(0).noalias() -= m1.row(0) * m1.adjoint();
  m2.row(0).noalias() -= m1.col(0).adjoint() * m1;
  m2.row(0).noalias() -= m1.col(0).adjoint() * m1.adjoint();
  VERIFY_IS_APPROX(m2,m2);
  
  m2.col(0).noalias() = m1.template triangularView<Upper>() * m1.col(0);
  m2.col(0).noalias() -= m1.adjoint().template triangularView<Upper>() * m1.col(0);
  m2.col(0).noalias() -= m1.template triangularView<Upper>() * m1.row(0).adjoint();
  m2.col(0).noalias() -= m1.adjoint().template triangularView<Upper>() * m1.row(0).adjoint();

  m2.row(0).noalias() = m1.row(0) * m1.template triangularView<Upper>();
  m2.row(0).noalias() -= m1.row(0) * m1.adjoint().template triangularView<Upper>();
  m2.row(0).noalias() -= m1.col(0).adjoint() * m1.template triangularView<Upper>();
  m2.row(0).noalias() -= m1.col(0).adjoint() * m1.adjoint().template triangularView<Upper>();
  VERIFY_IS_APPROX(m2,m2);
  
  m2.col(0).noalias() = m1.template selfadjointView<Upper>() * m1.col(0);
  m2.col(0).noalias() -= m1.adjoint().template selfadjointView<Upper>() * m1.col(0);
  m2.col(0).noalias() -= m1.template selfadjointView<Upper>() * m1.row(0).adjoint();
  m2.col(0).noalias() -= m1.adjoint().template selfadjointView<Upper>() * m1.row(0).adjoint();

  m2.row(0).noalias() = m1.row(0) * m1.template selfadjointView<Upper>();
  m2.row(0).noalias() -= m1.row(0) * m1.adjoint().template selfadjointView<Upper>();
  m2.row(0).noalias() -= m1.col(0).adjoint() * m1.template selfadjointView<Upper>();
  m2.row(0).noalias() -= m1.col(0).adjoint() * m1.adjoint().template selfadjointView<Upper>();
  VERIFY_IS_APPROX(m2,m2);
  
  m2.template selfadjointView<Lower>().rankUpdate(m1.col(0),-1);
  m2.template selfadjointView<Lower>().rankUpdate(m1.row(0),-1);

  // The following fancy matrix-matrix products are not safe yet regarding static allocation
//   m1 += m1.template triangularView<Upper>() * m2.col(;
//   m1.template selfadjointView<Lower>().rankUpdate(m2);
//   m1 += m1.template triangularView<Upper>() * m2;
//   m1 += m1.template selfadjointView<Lower>() * m2;
//   VERIFY_IS_APPROX(m1,m1);
}

template<typename Scalar>
void ctms_decompositions()
{
  const int maxSize = 16;
  const int size    = 12;

  typedef Eigen::Matrix<Scalar,
                        Eigen::Dynamic, Eigen::Dynamic,
                        0,
                        maxSize, maxSize> Matrix;

  typedef Eigen::Matrix<Scalar,
                        Eigen::Dynamic, 1,
                        0,
                        maxSize, 1> Vector;

  typedef Eigen::Matrix<std::complex<Scalar>,
                        Eigen::Dynamic, Eigen::Dynamic,
                        0,
                        maxSize, maxSize> ComplexMatrix;

  const Matrix A(Matrix::Random(size, size));
  const ComplexMatrix complexA(ComplexMatrix::Random(size, size));
  const Matrix saA = A.adjoint() * A;

  // Cholesky module
  Eigen::LLT<Matrix>  LLT;  LLT.compute(A);
  Eigen::LDLT<Matrix> LDLT; LDLT.compute(A);

  // Eigenvalues module
  Eigen::HessenbergDecomposition<ComplexMatrix> hessDecomp;        hessDecomp.compute(complexA);
  Eigen::ComplexSchur<ComplexMatrix>            cSchur(size);      cSchur.compute(complexA);
  Eigen::ComplexEigenSolver<ComplexMatrix>      cEigSolver;        cEigSolver.compute(complexA);
  Eigen::EigenSolver<Matrix>                    eigSolver;         eigSolver.compute(A);
  Eigen::SelfAdjointEigenSolver<Matrix>         saEigSolver(size); saEigSolver.compute(saA);
  Eigen::Tridiagonalization<Matrix>             tridiag;           tridiag.compute(saA);

  // LU module
  Eigen::PartialPivLU<Matrix> ppLU; ppLU.compute(A);
  Eigen::FullPivLU<Matrix>    fpLU; fpLU.compute(A);

  // QR module
  Eigen::HouseholderQR<Matrix>        hQR;  hQR.compute(A);
  Eigen::ColPivHouseholderQR<Matrix>  cpQR; cpQR.compute(A);
  Eigen::FullPivHouseholderQR<Matrix> fpQR; fpQR.compute(A);

  // SVD module
  Eigen::JacobiSVD<Matrix> jSVD; jSVD.compute(A, ComputeFullU | ComputeFullV);
}

void test_nomalloc()
{
  // check that our operator new is indeed called:
  VERIFY_RAISES_ASSERT(MatrixXd dummy(MatrixXd::Random(3,3)));
  CALL_SUBTEST_1(nomalloc(Matrix<float, 1, 1>()) );
  CALL_SUBTEST_2(nomalloc(Matrix4d()) );
  CALL_SUBTEST_3(nomalloc(Matrix<float,32,32>()) );
  
  // Check decomposition modules with dynamic matrices that have a known compile-time max size (ctms)
  CALL_SUBTEST_4(ctms_decompositions<float>());

}
