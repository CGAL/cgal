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
// Author(s)     : Mael Rouxel-Labbé
//                 Konstantinos Katrioplas
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
  typedef typename Traits::Vector                               Vector;

  CGAL_assertion(R.number_of_rows() == 3 && R.number_of_columns() == 3);
  CGAL_assertion(points.size() >= 3);

  FT xmin, ymin, zmin, xmax, ymax, zmax;
  xmin = ymin = zmin = FT{(std::numeric_limits<double>::max)()};
  xmax = ymax = zmax = FT{std::numeric_limits<double>::lowest()};

  for(const Point& pt : points)
  {
    Vector pv(3);
    pv.set(0, pt.x());
    pv.set(1, pt.y());
    pv.set(2, pt.z());
    pv = R * pv;

    xmin = (std::min)(xmin, pv(0));
    ymin = (std::min)(ymin, pv(1));
    zmin = (std::min)(zmin, pv(2));
    xmax = (std::max)(xmax, pv(0));
    ymax = (std::max)(ymax, pv(1));
    zmax = (std::max)(zmax, pv(2));
  }

  // volume
  return ((xmax - xmin) * (ymax - ymin) * (zmax - zmin));
}

// In this function, we only care about the fitness if it is smaller.
// Since it likely to not be smaller, we compute the volume and the value
// that is interesting from time to time.
// Since the volume only grows, we exit early if greater.
template <typename Traits, typename PointRange>
std::pair<typename Traits::FT, bool>
compute_fitness_if_smaller(const typename Traits::Matrix& R, // rotation matrix
                           const PointRange& points,
                           const typename Traits::FT cmp_val,
                           const Traits& /*traits*/)
{
  typedef typename Traits::FT                                   FT;
  typedef typename Traits::Point_3                              Point;
  typedef typename Traits::Vector                               Vector;

  CGAL_assertion(R.number_of_rows() == 3 && R.number_of_columns() == 3);
  CGAL_assertion(points.size() >= 3);

  FT xmin, ymin, zmin, xmax, ymax, zmax;
  xmin = ymin = zmin = FT{(std::numeric_limits<double>::max)()};
  xmax = ymax = zmax = FT{std::numeric_limits<double>::lowest()};

  // compute every 1%, with a minimum of 1000 iterations
  const std::size_t pn = points.size();
  const std::size_t mpn = (std::max)(pn / 10, std::size_t(1000));
  std::size_t next_check = mpn; // to avoid a costly modulo

//  std::cout << "pn, mpn, next_check: " << pn << " " << mpn << std::endl;

  std::size_t i = 0;
  for(const Point& pt : points)
  {
    Vector pv(3);
    pv.set(0, pt.x());
    pv.set(1, pt.y());
    pv.set(2, pt.z());
    pv = R * pv;

    xmin = (std::min)(xmin, pv(0));
    ymin = (std::min)(ymin, pv(1));
    zmin = (std::min)(zmin, pv(2));
    xmax = (std::max)(xmax, pv(0));
    ymax = (std::max)(ymax, pv(1));
    zmax = (std::max)(zmax, pv(2));

    if(i++ == next_check)
    {
      if((xmax - xmin) * (ymax - ymin) * (zmax - zmin) > cmp_val)
        return std::make_pair(FT(0), false);

      next_check += mpn;
    }
  }

  // volume
  return std::make_pair((xmax - xmin) * (ymax - ymin) * (zmax - zmin), true);
}

} // namespace internal
} // namespace Optimal_bounding_box
} // namespace CGAL

#endif // CGAL_OPTIMAL_BOUNDING_FITNESS_FUNCTION_H
