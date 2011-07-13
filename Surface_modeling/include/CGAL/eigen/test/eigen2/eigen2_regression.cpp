// This file is part of Eigen, a lightweight C++ template library
// for linear algebra. Eigen itself is part of the KDE project.
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
#include <Eigen/LeastSquares>

template<typename VectorType,
         typename HyperplaneType>
void makeNoisyCohyperplanarPoints(int numPoints,
                                  VectorType **points,
                                  HyperplaneType *hyperplane,
                                  typename VectorType::Scalar noiseAmplitude)
{
  typedef typename VectorType::Scalar Scalar;
  const int size = points[0]->size();
  // pick a random hyperplane, store the coefficients of its equation
  hyperplane->coeffs().resize(size + 1);
  for(int j = 0; j < size + 1; j++)
  {
    do {
      hyperplane->coeffs().coeffRef(j) = ei_random<Scalar>();
    } while(ei_abs(hyperplane->coeffs().coeff(j)) < 0.5);
  }

  // now pick numPoints random points on this hyperplane
  for(int i = 0; i < numPoints; i++)
  {
    VectorType& cur_point = *(points[i]);
    do
    {
      cur_point = VectorType::Random(size)/*.normalized()*/;
      // project cur_point onto the hyperplane
      Scalar x = - (hyperplane->coeffs().start(size).cwise()*cur_point).sum();
      cur_point *= hyperplane->coeffs().coeff(size) / x;
    } while( cur_point.norm() < 0.5
          || cur_point.norm() > 2.0 );
  }

  // add some noise to these points
  for(int i = 0; i < numPoints; i++ )
    *(points[i]) += noiseAmplitude * VectorType::Random(size);
}

template<typename VectorType>
void check_linearRegression(int numPoints,
                            VectorType **points,
                            const VectorType& original,
                            typename VectorType::Scalar tolerance)
{
  int size = points[0]->size();
  assert(size==2);
  VectorType result(size);
  linearRegression(numPoints, points, &result, 1);
  typename VectorType::Scalar error = (result - original).norm() / original.norm();
  VERIFY(ei_abs(error) < ei_abs(tolerance));
}

template<typename VectorType,
         typename HyperplaneType>
void check_fitHyperplane(int numPoints,
                         VectorType **points,
                         const HyperplaneType& original,
                         typename VectorType::Scalar tolerance)
{
  int size = points[0]->size();
  HyperplaneType result(size);
  fitHyperplane(numPoints, points, &result);
  result.coeffs() *= original.coeffs().coeff(size)/result.coeffs().coeff(size);
  typename VectorType::Scalar error = (result.coeffs() - original.coeffs()).norm() / original.coeffs().norm();
  std::cout << ei_abs(error) << "  xxx   " << ei_abs(tolerance) << std::endl;
  VERIFY(ei_abs(error) < ei_abs(tolerance));
}

void test_eigen2_regression()
{
  for(int i = 0; i < g_repeat; i++)
  {
#ifdef EIGEN_TEST_PART_1
    {
      Vector2f points2f [1000];
      Vector2f *points2f_ptrs [1000];
      for(int i = 0; i < 1000; i++) points2f_ptrs[i] = &(points2f[i]);
      Vector2f coeffs2f;
      Hyperplane<float,2> coeffs3f;
      makeNoisyCohyperplanarPoints(1000, points2f_ptrs, &coeffs3f, 0.01f);
      coeffs2f[0] = -coeffs3f.coeffs()[0]/coeffs3f.coeffs()[1];
      coeffs2f[1] = -coeffs3f.coeffs()[2]/coeffs3f.coeffs()[1];
      CALL_SUBTEST(check_linearRegression(10, points2f_ptrs, coeffs2f, 0.05f));
      CALL_SUBTEST(check_linearRegression(100, points2f_ptrs, coeffs2f, 0.01f));
      CALL_SUBTEST(check_linearRegression(1000, points2f_ptrs, coeffs2f, 0.002f));
    }
#endif
#ifdef EIGEN_TEST_PART_2
    {
      Vector2f points2f [1000];
      Vector2f *points2f_ptrs [1000];
      for(int i = 0; i < 1000; i++) points2f_ptrs[i] = &(points2f[i]);
      Hyperplane<float,2> coeffs3f;
      makeNoisyCohyperplanarPoints(1000, points2f_ptrs, &coeffs3f, 0.01f);
      CALL_SUBTEST(check_fitHyperplane(10, points2f_ptrs, coeffs3f, 0.05f));
      CALL_SUBTEST(check_fitHyperplane(100, points2f_ptrs, coeffs3f, 0.01f));
      CALL_SUBTEST(check_fitHyperplane(1000, points2f_ptrs, coeffs3f, 0.002f));
    }
#endif
#ifdef EIGEN_TEST_PART_3
    {
      Vector4d points4d [1000];
      Vector4d *points4d_ptrs [1000];
      for(int i = 0; i < 1000; i++) points4d_ptrs[i] = &(points4d[i]);
      Hyperplane<double,4> coeffs5d;
      makeNoisyCohyperplanarPoints(1000, points4d_ptrs, &coeffs5d, 0.01);
      CALL_SUBTEST(check_fitHyperplane(10, points4d_ptrs, coeffs5d, 0.05));
      CALL_SUBTEST(check_fitHyperplane(100, points4d_ptrs, coeffs5d, 0.01));
      CALL_SUBTEST(check_fitHyperplane(1000, points4d_ptrs, coeffs5d, 0.002));
    }
#endif
#ifdef EIGEN_TEST_PART_4
    {
      VectorXcd *points11cd_ptrs[1000];
      for(int i = 0; i < 1000; i++) points11cd_ptrs[i] = new VectorXcd(11);
      Hyperplane<std::complex<double>,Dynamic> *coeffs12cd = new Hyperplane<std::complex<double>,Dynamic>(11);
      makeNoisyCohyperplanarPoints(1000, points11cd_ptrs, coeffs12cd, 0.01);
      CALL_SUBTEST(check_fitHyperplane(100, points11cd_ptrs, *coeffs12cd, 0.025));
      CALL_SUBTEST(check_fitHyperplane(1000, points11cd_ptrs, *coeffs12cd, 0.006));
      delete coeffs12cd;
      for(int i = 0; i < 1000; i++) delete points11cd_ptrs[i];
    }
#endif
  }
}
