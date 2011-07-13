// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2008 Benoit Jacob <jacob.benoit.1@gmail.com>
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
#include <Eigen/StdVector>
#include <Eigen/Geometry>

template<typename MatrixType>
void check_stdvector_matrix(const MatrixType& m)
{
  typename MatrixType::Index rows = m.rows();
  typename MatrixType::Index cols = m.cols();
  MatrixType x = MatrixType::Random(rows,cols), y = MatrixType::Random(rows,cols);
  std::vector<MatrixType,Eigen::aligned_allocator<MatrixType> > v(10, MatrixType(rows,cols)), w(20, y);
  v[5] = x;
  w[6] = v[5];
  VERIFY_IS_APPROX(w[6], v[5]);
  v = w;
  for(int i = 0; i < 20; i++)
  {
    VERIFY_IS_APPROX(w[i], v[i]);
  }

  v.resize(21);
  v[20] = x;
  VERIFY_IS_APPROX(v[20], x);
  v.resize(22,y);
  VERIFY_IS_APPROX(v[21], y);
  v.push_back(x);
  VERIFY_IS_APPROX(v[22], x);
  VERIFY((size_t)&(v[22]) == (size_t)&(v[21]) + sizeof(MatrixType));

  // do a lot of push_back such that the vector gets internally resized
  // (with memory reallocation)
  MatrixType* ref = &w[0];
  for(int i=0; i<30 || ((ref==&w[0]) && i<300); ++i)
    v.push_back(w[i%w.size()]);
  for(unsigned int i=23; i<v.size(); ++i)
  {
    VERIFY(v[i]==w[(i-23)%w.size()]);
  }
}

template<typename TransformType>
void check_stdvector_transform(const TransformType&)
{
  typedef typename TransformType::MatrixType MatrixType;
  TransformType x(MatrixType::Random()), y(MatrixType::Random());
  std::vector<TransformType,Eigen::aligned_allocator<TransformType> > v(10), w(20, y);
  v[5] = x;
  w[6] = v[5];
  VERIFY_IS_APPROX(w[6], v[5]);
  v = w;
  for(int i = 0; i < 20; i++)
  {
    VERIFY_IS_APPROX(w[i], v[i]);
  }

  v.resize(21);
  v[20] = x;
  VERIFY_IS_APPROX(v[20], x);
  v.resize(22,y);
  VERIFY_IS_APPROX(v[21], y);
  v.push_back(x);
  VERIFY_IS_APPROX(v[22], x);
  VERIFY((size_t)&(v[22]) == (size_t)&(v[21]) + sizeof(TransformType));

  // do a lot of push_back such that the vector gets internally resized
  // (with memory reallocation)
  TransformType* ref = &w[0];
  for(int i=0; i<30 || ((ref==&w[0]) && i<300); ++i)
    v.push_back(w[i%w.size()]);
  for(unsigned int i=23; i<v.size(); ++i)
  {
    VERIFY(v[i].matrix()==w[(i-23)%w.size()].matrix());
  }
}

template<typename QuaternionType>
void check_stdvector_quaternion(const QuaternionType&)
{
  typedef typename QuaternionType::Coefficients Coefficients;
  QuaternionType x(Coefficients::Random()), y(Coefficients::Random());
  std::vector<QuaternionType,Eigen::aligned_allocator<QuaternionType> > v(10), w(20, y);
  v[5] = x;
  w[6] = v[5];
  VERIFY_IS_APPROX(w[6], v[5]);
  v = w;
  for(int i = 0; i < 20; i++)
  {
    VERIFY_IS_APPROX(w[i], v[i]);
  }

  v.resize(21);
  v[20] = x;
  VERIFY_IS_APPROX(v[20], x);
  v.resize(22,y);
  VERIFY_IS_APPROX(v[21], y);
  v.push_back(x);
  VERIFY_IS_APPROX(v[22], x);
  VERIFY((size_t)&(v[22]) == (size_t)&(v[21]) + sizeof(QuaternionType));

  // do a lot of push_back such that the vector gets internally resized
  // (with memory reallocation)
  QuaternionType* ref = &w[0];
  for(int i=0; i<30 || ((ref==&w[0]) && i<300); ++i)
    v.push_back(w[i%w.size()]);
  for(unsigned int i=23; i<v.size(); ++i)
  {
    VERIFY(v[i].coeffs()==w[(i-23)%w.size()].coeffs());
  }
}

void test_stdvector()
{
  // some non vectorizable fixed sizes
  CALL_SUBTEST_1(check_stdvector_matrix(Vector2f()));
  CALL_SUBTEST_1(check_stdvector_matrix(Matrix3f()));
  CALL_SUBTEST_2(check_stdvector_matrix(Matrix3d()));

  // some vectorizable fixed sizes
  CALL_SUBTEST_1(check_stdvector_matrix(Matrix2f()));
  CALL_SUBTEST_1(check_stdvector_matrix(Vector4f()));
  CALL_SUBTEST_1(check_stdvector_matrix(Matrix4f()));
  CALL_SUBTEST_2(check_stdvector_matrix(Matrix4d()));

  // some dynamic sizes
  CALL_SUBTEST_3(check_stdvector_matrix(MatrixXd(1,1)));
  CALL_SUBTEST_3(check_stdvector_matrix(VectorXd(20)));
  CALL_SUBTEST_3(check_stdvector_matrix(RowVectorXf(20)));
  CALL_SUBTEST_3(check_stdvector_matrix(MatrixXcf(10,10)));

  // some Transform
  CALL_SUBTEST_4(check_stdvector_transform(Projective2f()));
  CALL_SUBTEST_4(check_stdvector_transform(Projective3f()));
  CALL_SUBTEST_4(check_stdvector_transform(Projective3d()));
  //CALL_SUBTEST(heck_stdvector_transform(Projective4d()));

  // some Quaternion
  CALL_SUBTEST_5(check_stdvector_quaternion(Quaternionf()));
  CALL_SUBTEST_5(check_stdvector_quaternion(Quaterniond()));
}
