// Copyright (c) 2024 Max-Planck-Institute Saarbruecken (Germany), GeometryFactory (France)
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : André Nusser <anusser@mpi-inf.mpg.de>
//                 Marvin Künnemann <marvin@mpi-inf.mpg.de>
//                 Karl Bringmann <kbringma@mpi-inf.mpg.de>
//                 Andreas Fabri
// =============================================================================

#ifndef CGAL_FRECHET_DISTANCE_H
#define CGAL_FRECHET_DISTANCE_H

#include <CGAL/license/Frechet_distance.h>
#include <CGAL/basic.h>
#include <CGAL/Frechet_distance/internal/Frechet_distance.h>

#include <iterator>

namespace CGAL
{

/**
 * \ingroup PkgFrechetDistanceFunctions
 * decides if the Frechet distance between two polylines is larger than a given distance.
 *
 * \param polyline1 the first polyline defined by the sequence of consecutive points
 * \param polyline2 the second polyline defined by the sequence of consecutive points
 * \param distance the decision distance
 * \param traits the geometric traits object
 *
 * \tparam Traits a model of `FrechetDistanceTraits`
 * \tparam PointRange  a model of the concept `RandomAccessContainer`
 * with `Traits::Point` as value type.
 */
template < class Traits, class PointRange>
bool is_Frechet_distance_larger(const PointRange& polyline1,
                              const PointRange& polyline2,
                              const double distance,
                              const Traits& traits = Traits())
{
    Protect_FPU_rounding<true> p; // TODO only if using filtering!
    auto icurve1 = Frechet_distance_::internal::toCurve(polyline1, traits);
    auto icurve2 = Frechet_distance_::internal::toCurve(polyline2, traits);

    using distance_t = const typename decltype(icurve1)::distance_t;

    return ! Frechet_distance_::internal::lessThan(icurve1, icurve2, distance_t(distance));
}

/**
 * \ingroup PkgFrechetDistanceFunctions
 * approximates the Fréchet distance between two polylines up to an additive error
 * of `precision`.
 *
 * \param polyline1 the first polyline defined by the sequence of consecutive points
 * \param polyline2 the second polyline defined by the sequence of consecutive points
 * \param precision the precision of the approximation
 *  \param traits the geometric traits object
 *
 * \tparam Traits a model of `FrechetDistanceTraits`
 * \tparam PointRange  a model of the concept `RandomAccessContainer`
 * with `Traits::Point` as value type.
 *
 * @return an interval enclosing the exact result, the difference between the upper and
 * the lower bound being less than `precision`.
 */
template <class Traits,class PointRange>
std::pair<double,double> approximate_Frechet_distance(const PointRange& polyline1,
                                                      const PointRange& polyline2,
                                                      const double precision,
                                                      const Traits& traits = Traits())
{
    Protect_FPU_rounding<true> p; // TODO only if using filtering!
    auto icurve1 = Frechet_distance_::internal::toCurve(polyline1, traits);
    auto icurve2 = Frechet_distance_::internal::toCurve(polyline2, traits);

    return Frechet_distance_::internal::calcDistance(icurve1, icurve2, precision);
}

}  // end of namespace CGAL

#endif  // CGAL_FRECHET_DISTANCE_H
