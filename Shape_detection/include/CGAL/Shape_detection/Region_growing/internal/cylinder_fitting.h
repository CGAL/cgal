// Copyright (c) 2022 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Sven Oesau

// Header of the original source
// This file incorporates work covered by the following copyright and
// permission notice:
//
// David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2022
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt
// https://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
//
// Permission is hereby granted, free of charge, to any person or organization
// obtaining a copy of the software and accompanying documentation covered by
// this license(the "Software") to use, reproduce, display, distribute,
// execute, and transmit the Software, and to prepare derivative works of the
// Software, and to permit third - parties to whom the Software is furnished to
// do so, all subject to the following :
//
// The copyright notices in the Software and this entire statement, including
// the above license grant, this restriction and the following disclaimer,
// must be included in all copies of the Software, in whole or in part, and
// all derivative works of the Software, unless such copies or derivative
// works are solely in the form of machine - executable object code generated
// by a source language processor.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON - INFRINGEMENT.IN NO EVENT
// SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE
// FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE,
// ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.
//
// The code below uses the version of
// https://github.com/davideberly/GeometricTools available on 30th of May 2022.
// Version: 6.0.2022.01.06
//


#ifndef CGAL_SHAPE_DETECTION_INTERNAL_CYLINDER_FITTING_H
#define CGAL_SHAPE_DETECTION_INTERNAL_CYLINDER_FITTING_H

#include <cstdint>
#include <vector>
#include <thread>
#include <Eigen/Dense>
#include <CGAL/Origin.h>
#include <CGAL/license/Shape_detection.h>

// The algorithm for least-squares fitting of a point set by a cylinder is
// described in
//   https://www.geometrictools.com/Documentation/CylinderFitting.pdf

namespace CGAL
{
  namespace internal {
    template <class Kernel,
      typename Region,
      typename PointMap>
    void preprocess(const Region& indices,
      const PointMap point_map,
      typename Kernel::Vector_3& average,
      Eigen::VectorXd& mu,
      Eigen::MatrixXd& f0,
      Eigen::MatrixXd& f1,
      Eigen::MatrixXd& f2) {

      average = typename Kernel::Vector_3(0, 0, 0);

      for (auto item : indices) {
        average += get(point_map, item) - CGAL::ORIGIN;
      }

      average = average / static_cast<double>(indices.size());

      Eigen::MatrixXd prod(indices.size(), 6);
      prod.setZero();

      for (std::size_t i = 0; i < indices.size(); i++) {
        typename Kernel::Vector_3 p = get(point_map, indices[i])
                             - average - CGAL::ORIGIN;
        prod(i, 0) = p.x() * p.x();
        prod(i, 1) = p.x() * p.y();
        prod(i, 2) = p.x() * p.z();
        prod(i, 3) = p.y() * p.y();
        prod(i, 4) = p.y() * p.z();
        prod(i, 5) = p.z() * p.z();
        mu[0] += prod(i, 0);
        mu[1] += 2 * prod(i, 1);
        mu[2] += 2 * prod(i, 2);
        mu[3] += prod(i, 3);
        mu[4] += 2 * prod(i, 4);
        mu[5] += prod(i, 5);
      }

      mu /= static_cast<double>(indices.size());

      for (size_t i = 0; i < indices.size(); i++) {
        Eigen::VectorXd delta(6);
        delta[0] = prod(i, 0) - mu[0];
        delta[1] = 2.0 * prod(i, 1) - mu[1];
        delta[2] = 2.0 * prod(i, 2) - mu[2];
        delta[3] = prod(i, 3) - mu[3];
        delta[4] = 2.0 * prod(i, 4) - mu[4];
        delta[5] = prod(i, 5) - mu[5];

        f0(0, 0) += prod(i, 0);
        f0(0, 1) += prod(i, 1);
        f0(0, 2) += prod(i, 2);
        f0(1, 1) += prod(i, 3);
        f0(1, 2) += prod(i, 4);
        f0(2, 2) += prod(i, 5);

        typename Kernel::Vector_3 tmp = get(point_map, indices[i])
                               - average - CGAL::ORIGIN;
        Eigen::Vector3d v(tmp[0], tmp[1], tmp[2]);

        f1 += v * delta.transpose();
        f2 += delta * delta.transpose();
      }

      f0 = f0 * (1.0 / indices.size());
      f0(1, 0) = f0(0, 1);
      f0(2, 0) = f0(0, 2);
      f0(2, 1) = f0(1, 2);
      f1 = f1 * (1.0 / indices.size());
      f2 = f2 * (1.0 / indices.size());
    }

    template <class Kernel>
    typename Kernel::FT fit_cylinder(
      typename Kernel::Line_3& line,
      typename Kernel::FT& squared_radius,
      typename Kernel::Vector_3& direction,
      typename Kernel::Vector_3& average,
      typename Eigen::VectorXd& mu,
      typename Eigen::MatrixXd& f0,
      typename Eigen::MatrixXd& f1,
      typename Eigen::MatrixXd& f2) {

      Eigen::VectorXd dir(3);
      dir[0] = direction.x();
      dir[1] = direction.y();
      dir[2] = direction.z();

      Eigen::MatrixXd P(3, 3);
      P.setIdentity();
      P = P - (dir * dir.transpose());

      Eigen::MatrixXd S(3, 3);
      S(0, 0) = S(1, 1) = S(2, 2) = 0;
      S(0, 1) = -direction.z();
      S(0, 2) = direction.y();
      S(1, 0) = direction.z();
      S(1, 2) = -direction.x();
      S(2, 0) = -direction.y();
      S(2, 1) = direction.x();

      Eigen::MatrixXd A, Ahat, A2hat, Q;
      A = P * f0 * P;
      Ahat = -(S * A * S);
      A2hat = Ahat * A;

      double trace = A2hat(0, 0) + A2hat(1, 1) + A2hat(2, 2);
      Q = Ahat * (1.0 / trace);

      Eigen::VectorXd p(6), a, b;
      p[0] = P(0, 0);
      p[1] = P(0, 1);
      p[2] = P(0, 2);
      p[3] = P(1, 1);
      p[4] = P(1, 2);
      p[5] = P(2, 2);

      a = f1 * p;
      b = Q * a;

      double t0 = p.transpose() * (f2 * p);
      double t1 = 4.0 * a.transpose() * b;
      double t2 = 4.0 * b.transpose() * f0 * b;

      line = typename Kernel::Line_3(typename Kernel::Point_3(b[0] + average[0],
        b[1] + average[1], b[2] + average[2]), direction);

      squared_radius = (p.transpose() * mu) + double(b.transpose() * b);

      return t0 - t1 + t2;
    }
  }

  template<
    typename Kernel,
    typename Region,
    typename PointMap,
    typename NormalMap>
  typename Kernel::FT fit_cylinder(
    const Region& region, const PointMap point_map,
    const NormalMap normal_map,
    typename Kernel::Line_3& line,
    typename Kernel::FT& squared_radius,
    const Kernel &traits) {

    using FT = typename Kernel::FT;
    using Vector_3 = typename Kernel::Vector_3;
    typename Kernel::Construct_cross_product_vector_3 cross_product =
      traits.construct_cross_product_vector_3_object();

    squared_radius = -1.0;

    std::size_t nb = 0;
    Vector_3 mean_axis = CGAL::NULL_VECTOR;

    // Axis direction estimation from sample normals by averaging.
    for (std::size_t i = 0; i < region.size() - 1; ++i) {
      Vector_3 v0 = get(normal_map, region[i]);
      v0 = v0 / sqrt(v0 * v0);
      Vector_3 v1 = get(normal_map, region[i + 1]);
      v1 = v1 / sqrt(v1 * v1);
      Vector_3 axis = cross_product(v0, v1);
      if (sqrt(axis.squared_length()) < FT(1) / FT(100)) {
        continue;
      }
      axis = axis / sqrt(axis * axis);

      if (nb != 0 && (axis * mean_axis < 0)) {
        axis = -axis;
      }
      mean_axis = mean_axis + axis;
      ++nb;
    }

    // If no proper mean axis can be derived the fitting fails.
    if (mean_axis * mean_axis < FT(1) / FT(100))
      return (std::numeric_limits<double>::max)();

    mean_axis = mean_axis / sqrt(mean_axis * mean_axis);

    typename Kernel::Vector_3 average;

    Eigen::VectorXd mu(6);
    mu.setZero();
    Eigen::MatrixXd f0(3, 3);
    f0.setZero();
    Eigen::MatrixXd f1(3, 6);
    f1.setZero();
    Eigen::MatrixXd f2(6, 6);
    f2.setZero();

    internal::preprocess<Kernel, Region, PointMap>
      (region, point_map, average, mu, f0, f1, f2);

    return internal::fit_cylinder<Kernel>
      (line, squared_radius, mean_axis, average, mu, f0, f1, f2);
  }
}

#endif
