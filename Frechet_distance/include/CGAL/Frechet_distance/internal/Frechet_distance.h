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
#include <CGAL/Frechet_distance_traits_2.h>

namespace CGAL {
namespace Frechet_distance_ {
namespace internal {

template <class PointRange, class Traits>
auto toCurve(const PointRange& point_range, const Traits& /* traits */)
{
  using IPoint = typename Traits::Point;
  constexpr bool is_from_cgal_kernel =
    !std::is_same_v<typename Kernel_traits<IPoint>::Kernel, internal_kernel_traits::Dummy_kernel<IPoint>>;

  //TODO the traits should also be returned

  if constexpr (is_from_cgal_kernel)
  {
    if constexpr (std::is_floating_point<typename Traits::FT>::type::value &&
                  Kernel_traits<IPoint>::Kernel::Has_filtered_predicates_tag::value)
    {
      using AK = CGAL::Simple_cartesian<Interval_nt_advanced>;
      using EK = CGAL::Simple_cartesian<Exact_rational>;

      using Filtered_traits = std::pair<Frechet_distance_traits_2<AK>, Frechet_distance_traits_2<EK>>;

      return Curve<Filtered_traits, true>(point_range);
    }
  }
  return Curve<Traits, false>(point_range);
}

template <class Traits, bool is_filtered>
bool lessThan(Curve<Traits, is_filtered> const& curve1, Curve<Traits, is_filtered> const& curve2,
              const typename Curve<Traits, is_filtered>::distance_t& distance /*, const Traits& traits */)
{
  FrechetLight<Curve<Traits, is_filtered>> frechet;
  return frechet.lessThanWithFilters(distance, curve1, curve2);
}

template <typename Traits, bool is_filtered>
std::pair<double,double> calcDistance(Curve<Traits, is_filtered> const& curve1,
                                      Curve<Traits, is_filtered> const& curve2,
                                      double precision)
{
    FrechetLight<Curve<Traits, is_filtered>> frechet;
    return frechet.calcDistance(curve1, curve2, precision);
}

} } }  // end of namespace CGAL::Frechet_distance_::internal

#endif  // CGAL_INTERNAL_Frechet_distance_FRECHET_DISTANCE_H
