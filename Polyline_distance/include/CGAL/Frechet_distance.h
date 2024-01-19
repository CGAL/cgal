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
//
// =============================================================================

#ifndef CGAL_FRECHET_DISTANCE_H
#define CGAL_FRECHET_DISTANCE_H

#include <CGAL/Polyline_traits_2.h>
#include <CGAL/internal/Polyline_distance/Frechet_distance.h>

#include <CGAL/basic.h>

#include <iterator>

namespace CGAL {

/**
 * \ingroup PkgPolylineDistanceFunctions
 * decides if the Frechet distance between two polylines given as a range of points
 * is less than a given distance.
 * \param curve1 the first curve defined by the sequence of consecutive points along the polyline
 * \param curve2 the second curve defined by the sequence of consecutive points along the polyline
 * \param distance the decision distance
 * \tparam PointRange  a model of the concept `RandomAccessContainer`
 * with `Traits::Point` as value type.
 */
template <class PointRange,
          class Traits>
bool
continuous_Frechet_distance_less_than(const PointRange& curve1,
                                      const PointRange& curve2,
                                      const typename Traits::FT distance)
{
  auto icurve1 = toCurve(curve1);
  auto icurve2 = toCurve(curve2);
  auto idistance = toDistance(distance);

  return lessThan(icurve1, icurve2, idistance);
}

/**
 * \ingroup PkgPolylineDistanceFunctions
 * computes the Frechet distance between two polylines given as a range of points.
 * \param curve1 the first curve defined by the sequence of consecutive points along the polyline
 * \param curve2 the second curve defined by the sequence of consecutive points along the polyline
 * \tparam PointRange  a model of the concept `RandomAccessContainer`
 * with `Traits::Point` as value type.
 */
template <class PointRange,
          class Traits>
typename Traits::FT
continuous_Frechet_distance(const PointRange& curve1,
                            const PointRange& curve2,
							const double precision)
{
  auto icurve1 = toCurve(curve1);
  auto icurve2 = toCurve(curve2);

  return calcDistance(icurve1, icurve2);
}

} // end of namespace CGAL

#endif // CGAL_FRECHET_DISTANCE_H
