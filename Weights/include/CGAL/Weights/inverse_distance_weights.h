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

#ifndef CGAL_INVERSE_DISTANCE_WEIGHTS_H
#define CGAL_INVERSE_DISTANCE_WEIGHTS_H

#include <CGAL/Weights/internal/utils.h>

#include <CGAL/Point_2.h>
#include <CGAL/Point_3.h>

namespace CGAL {
namespace Weights {

/// \cond SKIP_IN_MANUAL
namespace inverse_distance_ns {

template<typename FT>
FT weight(const FT d)
{
  FT w = FT(0);
  CGAL_precondition(d != FT(0));
  if (d != FT(0))
    w = FT(1) / d;

  return w;
}

} // namespace inverse_distance_ns

template<typename GeomTraits>
typename GeomTraits::FT inverse_distance_weight(const typename GeomTraits::Point_2&,
                                                const typename GeomTraits::Point_2& r,
                                                const typename GeomTraits::Point_2&,
                                                const typename GeomTraits::Point_2& q,
                                                const GeomTraits& traits)
{
  using FT = typename GeomTraits::FT;

  const FT d = internal::distance_2(traits, q, r);
  return inverse_distance_ns::weight(d);
}

template<typename GeomTraits>
typename GeomTraits::FT inverse_distance_weight(const CGAL::Point_2<GeomTraits>& t,
                                                const CGAL::Point_2<GeomTraits>& r,
                                                const CGAL::Point_2<GeomTraits>& p,
                                                const CGAL::Point_2<GeomTraits>& q)
{
  const GeomTraits traits;
  return inverse_distance_weight(t, r, p, q, traits);
}

template<typename GeomTraits>
typename GeomTraits::FT inverse_distance_weight(const typename GeomTraits::Point_2& p,
                                                const typename GeomTraits::Point_2& q,
                                                const GeomTraits& traits)
{
  typename GeomTraits::Point_2 stub;
  return inverse_distance_weight(stub, p, stub, q, traits);
}

template<typename GeomTraits>
typename GeomTraits::FT inverse_distance_weight(const CGAL::Point_2<GeomTraits>& p,
                                                const CGAL::Point_2<GeomTraits>& q)
{
  CGAL::Point_2<GeomTraits> stub;
  return inverse_distance_weight(stub, p, stub, q);
}

template<typename GeomTraits>
typename GeomTraits::FT inverse_distance_weight(const typename GeomTraits::Point_3&,
                                                const typename GeomTraits::Point_3& r,
                                                const typename GeomTraits::Point_3&,
                                                const typename GeomTraits::Point_3& q,
                                                const GeomTraits& traits)
{
  using FT = typename GeomTraits::FT;

  const FT d = internal::distance_3(traits, q, r);
  return inverse_distance_ns::weight(d);
}

template<typename GeomTraits>
typename GeomTraits::FT inverse_distance_weight(const CGAL::Point_3<GeomTraits>& t,
                                                const CGAL::Point_3<GeomTraits>& r,
                                                const CGAL::Point_3<GeomTraits>& p,
                                                const CGAL::Point_3<GeomTraits>& q)
{
  const GeomTraits traits;
  return inverse_distance_weight(t, r, p, q, traits);
}

template<typename GeomTraits>
typename GeomTraits::FT inverse_distance_weight(const typename GeomTraits::Point_3& p,
                                                const typename GeomTraits::Point_3& q,
                                                const GeomTraits& traits)
{
  typename GeomTraits::Point_3 stub;
  return inverse_distance_weight(stub, p, stub, q, traits);
}

template<typename GeomTraits>
typename GeomTraits::FT inverse_distance_weight(const CGAL::Point_3<GeomTraits>& p,
                                                const CGAL::Point_3<GeomTraits>& q)
{
  CGAL::Point_3<GeomTraits> stub;
  return inverse_distance_weight(stub, p, stub, q);
}

/// \endcond

} // namespace Weights
} // namespace CGAL

#endif // CGAL_INVERSE_DISTANCE_WEIGHTS_H
