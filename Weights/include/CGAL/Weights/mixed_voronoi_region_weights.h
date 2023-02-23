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

#ifndef CGAL_MIXED_VORONOI_REGION_WEIGHTS_H
#define CGAL_MIXED_VORONOI_REGION_WEIGHTS_H

#include <CGAL/Weights/internal/utils.h>

#include <CGAL/Point_2.h>
#include <CGAL/Point_3.h>

namespace CGAL {
namespace Weights {

// 2D ==============================================================================================

/*!
  \ingroup PkgWeightsRefMixedVoronoiRegionWeights
  \brief computes the area of the mixed Voronoi cell in 2D using the points `p`, `q`, and `r`.
  \tparam GeomTraits a model of `AnalyticWeightTraits_2`
*/
template<typename GeomTraits>
typename GeomTraits::FT mixed_voronoi_area(const typename GeomTraits::Point_2& p,
                                           const typename GeomTraits::Point_2& q,
                                           const typename GeomTraits::Point_2& r,
                                           const GeomTraits& traits)
{
  using FT = typename GeomTraits::FT;
  using Point_2 = typename GeomTraits::Point_2;

  auto angle_2 = traits.angle_2_object();
  auto midpoint_2 = traits.construct_midpoint_2_object();
  auto circumcenter_2 = traits.construct_circumcenter_2_object();

  Point_2 center;
  if (angle_2(p, q, r) != CGAL::OBTUSE &&
      angle_2(q, r, p) != CGAL::OBTUSE &&
      angle_2(r, p, q) != CGAL::OBTUSE)
    center = circumcenter_2(p, q, r);
  else
    center = midpoint_2(r, p);

  const Point_2 m1 = midpoint_2(q, r);
  const Point_2 m2 = midpoint_2(q, p);

  const FT A1 = internal::positive_area_2(q, m1, center, traits);
  const FT A2 = internal::positive_area_2(q, center, m2, traits);

  return A1 + A2;
}

/*!
  \ingroup PkgWeightsRefMixedVoronoiRegionWeights
  \brief computes the area of the mixed Voronoi cell in 2D using the points `p`, `q`, and `r`.
  \tparam Kernel a model of `Kernel`
*/
template<typename GeomTraits>
typename GeomTraits::FT mixed_voronoi_area(const CGAL::Point_2<GeomTraits>& p,
                                           const CGAL::Point_2<GeomTraits>& q,
                                           const CGAL::Point_2<GeomTraits>& r)
{
  const GeomTraits traits;
  return mixed_voronoi_area(p, q, r, traits);
}

// 3D ==============================================================================================

/*!
  \ingroup PkgWeightsRefMixedVoronoiRegionWeights
  \brief computes the area of the mixed Voronoi cell in 3D using the points `p`, `q`, and `r`.
  \tparam GeomTraits a model of `AnalyticWeightTraits_3`
*/
template<typename GeomTraits>
typename GeomTraits::FT mixed_voronoi_area(const typename GeomTraits::Point_3& p,
                                           const typename GeomTraits::Point_3& q,
                                           const typename GeomTraits::Point_3& r,
                                           const GeomTraits& traits)
{
  using FT = typename GeomTraits::FT;
  using Point_3 = typename GeomTraits::Point_3;

  auto angle_3 = traits.angle_3_object();
  auto midpoint_3 = traits.construct_midpoint_3_object();
  auto circumcenter_3 = traits.construct_circumcenter_3_object();

  Point_3 center;
  if (angle_3(p, q, r) != CGAL::OBTUSE &&
      angle_3(q, r, p) != CGAL::OBTUSE &&
      angle_3(r, p, q) != CGAL::OBTUSE)
    center = circumcenter_3(p, q, r);
  else
    center = midpoint_3(r, p);

  const Point_3 m1 = midpoint_3(q, r);
  const Point_3 m2 = midpoint_3(q, p);

  const FT A1 = internal::positive_area_3(q, m1, center, traits);
  const FT A2 = internal::positive_area_3(q, center, m2, traits);

  return A1 + A2;
}

/*!
  \ingroup PkgWeightsRefMixedVoronoiRegionWeights
  \brief computes the area of the mixed Voronoi cell in 3D using the points `p`, `q`, and `r`.
  \tparam Kernel a model of `Kernel`
*/
template<typename Kernel>
typename Kernel::FT mixed_voronoi_area(const CGAL::Point_3<Kernel>& p,
                                       const CGAL::Point_3<Kernel>& q,
                                       const CGAL::Point_3<Kernel>& r)
{
  const Kernel traits;
  return mixed_voronoi_area(p, q, r, traits);
}

} // namespace Weights
} // namespace CGAL

#endif // CGAL_MIXED_VORONOI_REGION_WEIGHTS_H
