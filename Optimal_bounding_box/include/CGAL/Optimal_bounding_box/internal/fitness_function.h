// Copyright (c) 2018-2019 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Konstantinos Katrioplas
//                 Mael Rouxel-Labbé
//
#ifndef CGAL_OPTIMAL_BOUNDING_FITNESS_FUNCTION_H
#define CGAL_OPTIMAL_BOUNDING_FITNESS_FUNCTION_H

#include <CGAL/license/Optimal_bounding_box.h>

#include <CGAL/assertions.h>

#include <algorithm>
#include <limits>

namespace CGAL {
namespace Optimal_bounding_box {
namespace internal {

template <typename Traits, typename PointRange>
typename Traits::FT
compute_fitness(const typename Traits::Matrix& R, // rotation matrix
                const PointRange& points,
                const Traits& /*traits*/)
{
  typedef typename Traits::FT                                   FT;
  typedef typename Traits::Point_3                              Point;

  CGAL_assertion(R.number_of_rows() == 3 && R.number_of_columns() == 3);
  CGAL_assertion(points.size() >= 3);

  FT xmin, ymin, zmin, xmax, ymax, zmax;
  xmin = ymin = zmin = FT{std::numeric_limits<double>::max()};
  xmax = ymax = zmax = FT{std::numeric_limits<double>::lowest()};

  for(const Point& pt : points)
  {
    const FT x = pt.x(), y = pt.y(), z = pt.z();

    const FT rx = x*R(0, 0) + y*R(0, 1) + z*R(0, 2);
    const FT ry = x*R(1, 0) + y*R(1, 1) + z*R(1, 2);
    const FT rz = x*R(2, 0) + y*R(2, 1) + z*R(2, 2);

    xmin = (std::min)(xmin, rx);
    ymin = (std::min)(ymin, ry);
    zmin = (std::min)(zmin, rz);
    xmax = (std::max)(xmax, rx);
    ymax = (std::max)(ymax, ry);
    zmax = (std::max)(zmax, rz);
  }

  // volume
  return ((xmax - xmin) * (ymax - ymin) * (zmax - zmin));
}

} // namespace internal
} // namespace Optimal_bounding_box
} // namespace CGAL

#endif // CGAL_OPTIMAL_BOUNDING_FITNESS_FUNCTION_H
