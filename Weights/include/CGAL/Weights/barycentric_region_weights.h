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

#ifndef CGAL_BARYCENTRIC_REGION_WEIGHTS_H
#define CGAL_BARYCENTRIC_REGION_WEIGHTS_H

#include <CGAL/Weights/utils.h>

#include <CGAL/Point_2.h>
#include <CGAL/Point_3.h>

namespace CGAL {
namespace Weights {

template<typename GeomTraits>
typename GeomTraits::FT barycentric_area(const typename GeomTraits::Point_2& p,
                                         const typename GeomTraits::Point_2& q,
                                         const typename GeomTraits::Point_2& r,
                                         const GeomTraits& traits)
{
  using FT = typename GeomTraits::FT;
  using Point_2 = typename GeomTraits::Point_2;

  auto midpoint_2 = traits.construct_midpoint_2_object();
  auto centroid_2 = traits.construct_centroid_2_object();

  const Point_2 center = centroid_2(p, q, r);
  const Point_2 m1 = midpoint_2(q, r);
  const Point_2 m2 = midpoint_2(q, p);

  const FT A1 = internal::positive_area_2(traits, q, m1, center);
  const FT A2 = internal::positive_area_2(traits, q, center, m2);
  return A1 + A2;
}

template<typename GeomTraits>
typename GeomTraits::FT barycentric_area(const CGAL::Point_2<GeomTraits>& p,
                                         const CGAL::Point_2<GeomTraits>& q,
                                         const CGAL::Point_2<GeomTraits>& r)
{
  const GeomTraits traits;
  return barycentric_area(p, q, r, traits);
}

template<typename GeomTraits>
typename GeomTraits::FT barycentric_area(const typename GeomTraits::Point_3& p,
                                         const typename GeomTraits::Point_3& q,
                                         const typename GeomTraits::Point_3& r,
                                         const GeomTraits& traits)
{
  using FT = typename GeomTraits::FT;
  using Point_3 = typename GeomTraits::Point_3;

  auto midpoint_3 = traits.construct_midpoint_3_object();
  auto centroid_3 = traits.construct_centroid_3_object();

  const Point_3 center = centroid_3(p, q, r);
  const Point_3 m1 = midpoint_3(q, r);
  const Point_3 m2 = midpoint_3(q, p);

  const FT A1 = internal::positive_area_3(traits, q, m1, center);
  const FT A2 = internal::positive_area_3(traits, q, center, m2);
  return A1 + A2;
}

template<typename GeomTraits>
typename GeomTraits::FT barycentric_area(const CGAL::Point_3<GeomTraits>& p,
                                         const CGAL::Point_3<GeomTraits>& q,
                                         const CGAL::Point_3<GeomTraits>& r)
{
  const GeomTraits traits;
  return barycentric_area(p, q, r, traits);
}

} // namespace Weights
} // namespace CGAL

#endif // CGAL_BARYCENTRIC_REGION_WEIGHTS_H
