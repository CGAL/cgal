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

#ifndef CGAL_SHEPARD_WEIGHTS_H
#define CGAL_SHEPARD_WEIGHTS_H

#include <CGAL/Weights/internal/utils.h>

#include <CGAL/Point_2.h>
#include <CGAL/Point_3.h>

namespace CGAL {
namespace Weights {

/// \cond SKIP_IN_MANUAL

namespace shepard_ns {

template<typename GeomTraits>
typename GeomTraits::FT weight(const GeomTraits& traits,
                               const typename GeomTraits::FT d,
                               const typename GeomTraits::FT p)
{
  using FT = typename GeomTraits::FT;

  FT w = FT(0);
  CGAL_precondition(is_positive(d));
  if(is_positive(d))
    w = internal::power(d, -p);

  return w;
}

} // namespace shepard_ns

template<typename GeomTraits>
typename GeomTraits::FT shepard_weight(const typename GeomTraits::Point_2&,
                                       const typename GeomTraits::Point_2& r,
                                       const typename GeomTraits::Point_2&,
                                       const typename GeomTraits::Point_2& q,
                                       const typename GeomTraits::FT a,
                                       const GeomTraits& traits)
{
  using FT = typename GeomTraits::FT;

  const FT d = internal::distance_2(traits, q, r);
  return shepard_ns::weight(traits, d, a);
}

template<typename GeomTraits>
typename GeomTraits::FT shepard_weight(const CGAL::Point_2<GeomTraits>& t,
                                       const CGAL::Point_2<GeomTraits>& r,
                                       const CGAL::Point_2<GeomTraits>& p,
                                       const CGAL::Point_2<GeomTraits>& q,
                                       const typename GeomTraits::FT a = typename GeomTraits::FT(1))
{
  const GeomTraits traits;
  return shepard_weight(t, r, p, q, a, traits);
}

template<typename GeomTraits>
typename GeomTraits::FT shepard_weight(const typename GeomTraits::Point_2& p,
                                       const typename GeomTraits::Point_2& q,
                                       const typename GeomTraits::FT a,
                                       const GeomTraits& traits)
{
  typename GeomTraits::Point_2 stub;
  return shepard_weight(stub, p, stub, q, a, traits);
}

template<typename GeomTraits>
typename GeomTraits::FT shepard_weight(const CGAL::Point_2<GeomTraits>& p,
                                       const CGAL::Point_2<GeomTraits>& q,
                                       const typename GeomTraits::FT a = typename GeomTraits::FT(1))
{
  CGAL::Point_2<GeomTraits> stub;
  return shepard_weight(stub, p, stub, q, a);
}

template<typename GeomTraits>
typename GeomTraits::FT shepard_weight(const typename GeomTraits::Point_3&,
                                       const typename GeomTraits::Point_3& r,
                                       const typename GeomTraits::Point_3&,
                                       const typename GeomTraits::Point_3& q,
                                       const typename GeomTraits::FT a,
                                       const GeomTraits& traits)
{
  using FT = typename GeomTraits::FT;
  const FT d = internal::distance_3(traits, q, r);
  return shepard_ns::weight(traits, d, a);
}

template<typename GeomTraits>
typename GeomTraits::FT shepard_weight(const CGAL::Point_3<GeomTraits>& t,
                                       const CGAL::Point_3<GeomTraits>& r,
                                       const CGAL::Point_3<GeomTraits>& p,
                                       const CGAL::Point_3<GeomTraits>& q,
                                       const typename GeomTraits::FT a = typename GeomTraits::FT(1))
{
  const GeomTraits traits;
  return shepard_weight(t, r, p, q, a, traits);
}

template<typename GeomTraits>
typename GeomTraits::FT shepard_weight(const typename GeomTraits::Point_3& p,
                                       const typename GeomTraits::Point_3& q,
                                       const typename GeomTraits::FT a,
                                       const GeomTraits& traits)
{
  typename GeomTraits::Point_3 stub;
  return shepard_weight(stub, p, stub, q, a, traits);
}

template<typename GeomTraits>
typename GeomTraits::FT shepard_weight(const CGAL::Point_3<GeomTraits>& p,
                                       const CGAL::Point_3<GeomTraits>& q,
                                       const typename GeomTraits::FT a = typename GeomTraits::FT(1))
{
  CGAL::Point_3<GeomTraits> stub;
  return shepard_weight(stub, p, stub, q, a);
}

/// \endcond

} // namespace Weights
} // namespace CGAL

#endif // CGAL_SHEPARD_WEIGHTS_H
