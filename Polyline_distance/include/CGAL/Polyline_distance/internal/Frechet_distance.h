// Copyright (c) 2019 Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
//
// Author(s)     : André Nusser <anusser@mpi-inf.mpg.de>
//                 Marvin Künnemann <marvin@mpi-inf.mpg.de>
//                 Karl Bringmann <kbringma@mpi-inf.mpg.de>
//
// =============================================================================

#ifndef CGAL_INTERNAL_POLYLINE_DISTANCE_FRECHET_DISTANCE_H
#define CGAL_INTERNAL_POLYLINE_DISTANCE_FRECHET_DISTANCE_H
#include <CGAL/license/Polyline_distance.h>
#include <CGAL/Polyline_distance/internal/curve.h>
#include <CGAL/Polyline_distance/internal/frechet_light.h>
#include <CGAL/Polyline_distance/internal/geometry_basics.h>

namespace CGAL {
namespace Polyline_distance {
namespace internal {

  template <class PointRange, class Traits>
  auto toCurve(const PointRange& point_range, const Traits& traits)
{
    typedef typename PointRange::const_iterator iterator;
    typedef typename std::iterator_traits<iterator>::value_type Point;
    typedef typename CGAL::Kernel_traits<Point>::Kernel K;

    Curve<Traits> curve(point_range);

    return curve;
}

template <class NT>
//distance_t 
auto toDistance(NT distance)
{
    return to_interval(distance);
}


template <class Traits>
bool lessThan(Curve<Traits> const& curve1, Curve<Traits> const& curve2,
              const typename Traits::distance_t& distance, const Traits& traits)
{
  FrechetLight<Traits> frechet;
  return frechet.lessThanWithFilters(distance, curve1, curve2);
}

template <typename Traits>
typename Traits::distance_t calcDistance(Curve<Traits> const& curve1, Curve<Traits> const& curve2)
{
    FrechetLight<Traits> frechet;
    return frechet.calcDistance(curve1, curve2);
}

} // namespace internal
} // namespace Polyline_distance
}  // end of namespace CGAL

#endif  // CGAL_INTERNAL_POLYLINE_DISTANCE_FRECHET_DISTANCE_H
