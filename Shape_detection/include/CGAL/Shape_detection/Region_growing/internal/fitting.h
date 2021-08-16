// Copyright (c) 2021 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Simon Giraudot
//

#ifndef CGAL_SHAPE_DETECTION_REGION_GROWING_INTERNAL_FITTING_H
#define CGAL_SHAPE_DETECTION_REGION_GROWING_INTERNAL_FITTING_H

#include <CGAL/license/Shape_detection.h>

#include <CGAL/Eigen_diagonalize_traits.h>
#include <CGAL/property_map.h>

namespace CGAL {
namespace Shape_detection {
namespace internal {

template <typename PointRange,
          typename Sqrt, typename Squared_distance_2,
          typename Point_2, typename FT>
bool circle_fit (const PointRange& points,
                 const Sqrt& sqrt,
                 const Squared_distance_2& squared_distance_2,
                 Point_2& center, FT& radius)
{
  using Diagonalize_traits = Eigen_diagonalize_traits<FT, 4>;
  using Covariance_matrix = typename Diagonalize_traits::Covariance_matrix;
  using Vector = typename Diagonalize_traits::Vector;
  using Matrix = typename Diagonalize_traits::Matrix;

  // Use bbox to compute diagonalization with smaller coordinates to
  // avoid loss of precision when inverting large coordinates
  Bbox_2 bbox = bbox_2 (points.begin(), points.end());

  // Circle least squares fitting,
  // Circle of center (a,b) and radius R
  // Ri = sqrt((xi - a)^2 + (yi - b)^2)
  // Minimize Sum(Ri^2 - R^2)^2
  // -> Minimize Sum(xi^2 + yi^2 − 2 a*xi − 2 b*yi + a^2 + b^2 − R^2)^2
  // let B=-2a ; C=-2b; D= a^2 + b^2 - R^2
  // let ri = x^2 + y^2
  // -> Minimize Sum(D + B*xi + C*yi + ri)^2
  // -> Minimize Sum(1 + B/D*xi + C/D*yi + ri/D)^2
  // -> system of linear equations
  // -> diagonalize matrix
  //    NB   x   y   r
  //        xx  xy  xr
  //            yy  yr
  //                rr
  //
  // -> center coordinates = -0.5 * eigenvector(1) / eigenvector(3) ; -0.5 * eigenvector(2) / eigenvector(3)
  Covariance_matrix A
    = { FT(0), FT(0), FT(0), FT(0), FT(0),
        FT(0), FT(0), FT(0), FT(0), FT(0) };

  A[0] = static_cast<FT>(points.size());
  for (const Point_2& p : points)
  {
    FT x = p.x() - bbox.xmin();
    FT y = p.y() - bbox.ymin();
    FT r = x*x + y*y;
    A[1] += x;
    A[2] += y;
    A[3] += r;
    A[4] += x * x;
    A[5] += x * y;
    A[6] += x * r;
    A[7] += y * y;
    A[8] += y * r;
    A[9] += r * r;
  }

  Vector eigenvalues = { FT(0), FT(0), FT(0), FT(0) };
  Matrix eigenvectors = { FT(0), FT(0), FT(0), FT(0),
                          FT(0), FT(0), FT(0), FT(0),
                          FT(0), FT(0), FT(0), FT(0),
                          FT(0), FT(0), FT(0), FT(0) };

  Diagonalize_traits::diagonalize_selfadjoint_covariance_matrix
    (A, eigenvalues, eigenvectors);

  // Perfect line case, no circle can be fitted
  if (eigenvectors[3] == 0)
    return false;

  center = Point_2 (bbox.xmin() - FT(0.5) * (eigenvectors[1] / eigenvectors[3]),
                    bbox.ymin() - FT(0.5) * (eigenvectors[2] / eigenvectors[3]));

  radius = FT(0);
  for (const Point_2& p : points)
    radius += sqrt (squared_distance_2 (p, center));
  radius /= points.size();

  return true;
}

template <typename PointRange,
          typename Sqrt, typename Squared_distance_3,
          typename Point_3, typename FT>
bool sphere_fit (const PointRange& points,
                 const Sqrt& sqrt,
                 const Squared_distance_3& squared_distance_3,
                 Point_3& center, FT& radius)
{
  using Diagonalize_traits = Eigen_diagonalize_traits<FT, 5>;
  using Covariance_matrix = typename Diagonalize_traits::Covariance_matrix;
  using Vector = typename Diagonalize_traits::Vector;
  using Matrix = typename Diagonalize_traits::Matrix;

  // Use bbox to compute diagonalization with smaller coordinates to
  // avoid loss of precision when inverting large coordinates
  Bbox_3 bbox = bbox_3 (points.begin(), points.end());

  // Sphere least squares fitting
  // (see Least_square_circle_fit_region for details about computation)
  Covariance_matrix A
    = { FT(0), FT(0), FT(0), FT(0), FT(0),
        FT(0), FT(0), FT(0), FT(0), FT(0),
        FT(0), FT(0), FT(0), FT(0), FT(0) };

  A[0] = static_cast<FT>(points.size());
  for (const Point_3& p : points)
  {
    FT x = p.x() - bbox.xmin();
    FT y = p.y() - bbox.ymin();
    FT z = p.z() - bbox.zmin();
    FT r = x*x + y*y + z*z;
    A[1] += x;
    A[2] += y;
    A[3] += z;
    A[4] += r;
    A[5] += x * x;
    A[6] += x * y;
    A[7] += x * z;
    A[8] += x * r;
    A[9] += y * y;
    A[10] += y * z;
    A[11] += y * r;
    A[12] += z * z;
    A[13] += z * r;
    A[14] += r * r;
  }

  Vector eigenvalues = { FT(0), FT(0), FT(0), FT(0), FT(0) };
  Matrix eigenvectors = { FT(0), FT(0), FT(0), FT(0), FT(0),
                          FT(0), FT(0), FT(0), FT(0), FT(0),
                          FT(0), FT(0), FT(0), FT(0), FT(0),
                          FT(0), FT(0), FT(0), FT(0), FT(0),
                          FT(0), FT(0), FT(0), FT(0), FT(0) };

  Diagonalize_traits::diagonalize_selfadjoint_covariance_matrix
    (A, eigenvalues, eigenvectors);

  // Perfect plane case, no sphere can be fitted
  if (eigenvectors[4] == 0)
    return false;

  center = Point_3 (bbox.xmin() - FT(0.5) * (eigenvectors[1] / eigenvectors[4]),
                    bbox.ymin() - FT(0.5) * (eigenvectors[2] / eigenvectors[4]),
                    bbox.zmin() - FT(0.5) * (eigenvectors[3] / eigenvectors[4]));

  radius = FT(0);
  for (const Point_3& p : points)
    radius += sqrt (squared_distance_3 (p, center));
  radius /= points.size();

  return true;
}

template <typename PointRange, typename PointMap, typename NormalMap,
          typename Sqrt, typename Squared_distance_3,
          typename Line_3, typename FT>
bool cylinder_fit (const PointRange& points,
                   PointMap point_map, NormalMap normal_map,
                   const Sqrt& sqrt,
                   const Squared_distance_3& squared_distance_3,
                   Line_3& axis, FT& radius)
{
  using Point_3 = typename boost::property_traits<PointMap>::value_type;
  using Vector_3 = typename boost::property_traits<NormalMap>::value_type;

  const Point_3& ref = get(point_map, *(points.begin()));

  Vector_3 mean_axis = CGAL::NULL_VECTOR;
  std::size_t nb = 0;
  radius = FT(0);
  Point_3 point_on_axis = ORIGIN;
  for (std::size_t i = 0; i < points.size() - 1; ++ i)
  {
    Vector_3 v0 = get(normal_map, *std::next(points.begin(), i));
    v0 = v0 / sqrt(v0*v0);
    Vector_3 v1 = get(normal_map, *std::next(points.begin(), i + 1));
    v1 = v1 / sqrt(v1*v1);
    Vector_3 axis = cross_product (v0, v1);
    if (sqrt(axis.squared_length()) < (FT)(0.01))
      continue;
    axis = axis / sqrt(axis * axis);

    const Point_3& p0 = get(point_map, *std::next(points.begin(), i));
    const Point_3& p1 = get(point_map, *std::next(points.begin(), i + 1));

    Vector_3 xdir = v0 - axis * (v0 * axis);
    xdir = xdir / sqrt (xdir * xdir);

    Vector_3 ydir = cross_product (axis, xdir);
    ydir = ydir / sqrt (ydir * ydir);

    FT v1x = v1 * ydir;
    FT v1y = -v1 * xdir;

    Vector_3 d (p0, p1);

    FT ox = xdir * d;
    FT oy = ydir * d;
    FT ldist = v1x * ox + v1y * oy;

    FT radius = ldist / v1x;

    Point_3 point = p0 + xdir * radius;
    Line_3 line (point, axis);
    point = line.projection(ref);

    point_on_axis = barycenter (point_on_axis, static_cast<FT>(nb), point, FT(1));

    radius += abs(radius);

    if (nb != 0 && axis * mean_axis < 0)
      axis = -axis;

    mean_axis = mean_axis + axis;
    ++ nb;
  }

  if (nb == 0)
    return false;

  mean_axis = mean_axis / sqrt(mean_axis * mean_axis);
  axis = Line_3 (point_on_axis, mean_axis);
  radius /= nb;

  radius = FT(0);
  for (const auto& p : points)
  {
    const Point_3& p0 = get(point_map, p);
    radius += sqrt(squared_distance_3(p0, axis));
  }
  radius /= points.size();

  return true;
}

} // internal
} // namespace Shape_detection
} // namespace CGAL


#endif // CGAL_SHAPE_DETECTION_REGION_GROWING_POINT_SET_INTERNAL_FITTING_H
