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

#ifndef CGAL_INTERNAL_Frechet_distance_FRECHET_DISTANCE_H
#define CGAL_INTERNAL_Frechet_distance_FRECHET_DISTANCE_H
#include <CGAL/license/Frechet_distance.h>
#include <CGAL/Frechet_distance/internal/curve.h>
#include <CGAL/Frechet_distance/internal/frechet_light.h>
#include <CGAL/Frechet_distance/internal/geometry_basics.h>

namespace CGAL {
namespace Frechet_distance_ {
namespace internal {

template <class PointRange, class Traits>
auto toCurve(const PointRange& point_range, const Traits& /* traits */)
{
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
              const typename Curve<Traits>::distance_t& distance, const Traits& /* traits */)
{
  FrechetLight<Curve<Traits>> frechet;
  return frechet.lessThanWithFilters(distance, curve1, curve2);
}

template <typename Traits>
std::pair<double,double> calcDistance(Curve<Traits> const& curve1,
                                      Curve<Traits> const& curve2,
                                      double precision)
{
    FrechetLight<Curve<Traits>> frechet;
    return frechet.calcDistance(curve1, curve2, precision);
}

} } }  // end of namespace CGAL::Frechet_distance_::internal

#endif  // CGAL_INTERNAL_Frechet_distance_FRECHET_DISTANCE_H
