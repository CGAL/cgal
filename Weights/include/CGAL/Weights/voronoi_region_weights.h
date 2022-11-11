// Copyright (c) 2020 GeometryFactory SARL (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Dmitry Anisimov
//

#ifndef CGAL_VORONOI_REGION_WEIGHTS_H
#define CGAL_VORONOI_REGION_WEIGHTS_H

// Internal includes.
#include <CGAL/Weights/internal/utils.h>

namespace CGAL {
namespace Weights {

  #if defined(DOXYGEN_RUNNING)

  /*!
    \ingroup PkgWeightsRefVoronoiRegionWeights

    \brief computes the area of the Voronoi cell in 2D using the points `p`, `q`
    and `r`, given a traits class `traits` with geometric objects, predicates, and constructions.
  */
  template<typename GeomTraits>
  typename GeomTraits::FT voronoi_area(
    const typename GeomTraits::Point_2& p,
    const typename GeomTraits::Point_2& q,
    const typename GeomTraits::Point_2& r,
    const GeomTraits& traits) { }

  /*!
    \ingroup PkgWeightsRefVoronoiRegionWeights

    \brief computes the area of the Voronoi cell in 3D using the points `p`, `q`
    and `r`, given a traits class `traits` with geometric objects, predicates, and constructions.
  */
  template<typename GeomTraits>
  typename GeomTraits::FT voronoi_area(
    const typename GeomTraits::Point_3& p,
    const typename GeomTraits::Point_3& q,
    const typename GeomTraits::Point_3& r,
    const GeomTraits& traits) { }

  /*!
    \ingroup PkgWeightsRefVoronoiRegionWeights

    \brief computes the area of the Voronoi cell in 2D using the points `p`, `q`
    and `r` which are parameterized by a `Kernel` K.
  */
  template<typename K>
  typename K::FT voronoi_area(
    const CGAL::Point_2<K>& p,
    const CGAL::Point_2<K>& q,
    const CGAL::Point_2<K>& r) { }

  /*!
    \ingroup PkgWeightsRefVoronoiRegionWeights

    \brief computes the area of the Voronoi cell in 3D using the points `p`, `q`
    and `r` which are parameterized by a `Kernel` K.
  */
  template<typename K>
  typename K::FT voronoi_area(
    const CGAL::Point_3<K>& p,
    const CGAL::Point_3<K>& q,
    const CGAL::Point_3<K>& r) { }

  #endif // DOXYGEN_RUNNING

  /// \cond SKIP_IN_MANUAL
  template<typename GeomTraits>
  typename GeomTraits::FT voronoi_area(
    const typename GeomTraits::Point_2& p,
    const typename GeomTraits::Point_2& q,
    const typename GeomTraits::Point_2& r,
    const GeomTraits& traits) {

    using FT = typename GeomTraits::FT;
    const auto circumcenter_2 =
      traits.construct_circumcenter_2_object();
    const auto midpoint_2 =
      traits.construct_midpoint_2_object();

    const auto center = circumcenter_2(p, q, r);
    const auto m1 = midpoint_2(q, r);
    const auto m2 = midpoint_2(q, p);

    const FT A1 = internal::positive_area_2(traits, q, m1, center);
    const FT A2 = internal::positive_area_2(traits, q, center, m2);
    return A1 + A2;
  }

  template<typename GeomTraits>
  typename GeomTraits::FT voronoi_area(
    const CGAL::Point_2<GeomTraits>& p,
    const CGAL::Point_2<GeomTraits>& q,
    const CGAL::Point_2<GeomTraits>& r) {

    const GeomTraits traits;
    return voronoi_area(p, q, r, traits);
  }

  template<typename GeomTraits>
  typename GeomTraits::FT voronoi_area(
    const typename GeomTraits::Point_3& p,
    const typename GeomTraits::Point_3& q,
    const typename GeomTraits::Point_3& r,
    const GeomTraits& traits) {

    using FT = typename GeomTraits::FT;
    const auto circumcenter_3 =
      traits.construct_circumcenter_3_object();
    const auto midpoint_3 =
      traits.construct_midpoint_3_object();

    const auto center = circumcenter_3(p, q, r);
    const auto m1 = midpoint_3(q, r);
    const auto m2 = midpoint_3(q, p);

    const FT A1 = internal::positive_area_3(traits, q, m1, center);
    const FT A2 = internal::positive_area_3(traits, q, center, m2);
    return A1 + A2;
  }

  template<typename GeomTraits>
  typename GeomTraits::FT voronoi_area(
    const CGAL::Point_3<GeomTraits>& p,
    const CGAL::Point_3<GeomTraits>& q,
    const CGAL::Point_3<GeomTraits>& r) {

    const GeomTraits traits;
    return voronoi_area(p, q, r, traits);
  }
  /// \endcond

} // namespace Weights
} // namespace CGAL

#endif // CGAL_VORONOI_REGION_WEIGHTS_H
