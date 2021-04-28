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

namespace CGAL {
namespace Shape_detection {
namespace internal {

template <typename PointRange,
          typename Sqrt, typename Squared_distance_3,
          typename Point_3, typename FT>
void sphere_fit (const PointRange& points,
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

  A[0] = points.size();
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

  center = Point_3 (bbox.xmin() - FT(0.5) * (eigenvectors[1] / eigenvectors[4]),
                    bbox.ymin() - FT(0.5) * (eigenvectors[2] / eigenvectors[4]),
                    bbox.zmin() - FT(0.5) * (eigenvectors[3] / eigenvectors[4]));

  radius = FT(0);
  for (const Point_3& p : points)
    radius += sqrt (squared_distance_3 (p, center));
  radius /= points.size();
}


} // internal
} // namespace Shape_detection
} // namespace CGAL


#endif // CGAL_SHAPE_DETECTION_REGION_GROWING_POINT_SET_INTERNAL_FITTING_H
