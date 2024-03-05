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

#include <CGAL/Polyline_distance/internal/curve.h>
#include <CGAL/Polyline_distance/internal/frechet_light.h>
#include <CGAL/Polyline_distance/internal/geometry_basics.h>

namespace CGAL {
namespace Polyline_distance {
namespace internal {

template <class PointRange>
Curve toCurve(const PointRange& point_range)
{
    Curve curve;

    curve.reserve(point_range.size());
    for (auto const& point : point_range) {
        auto ipoint = Point(point.x(), point.y());
        curve.push_back(ipoint);
    }

    return curve;
}

template <class NT>
distance_t toDistance(NT distance)
{
    return distance;
}

bool lessThan(Curve const& curve1, Curve const& curve2, distance_t distance)
{
    FrechetLight frechet;
    return frechet.lessThanWithFilters(distance, curve1, curve2);
}

distance_t calcDistance(Curve const& curve1, Curve const& curve2)
{
    FrechetLight frechet;
    return frechet.calcDistance(curve1, curve2);
}

} // namespace internal
} // namespace Polyline_distance
}  // end of namespace CGAL

#endif  // CGAL_INTERNAL_POLYLINE_DISTANCE_FRECHET_DISTANCE_H
