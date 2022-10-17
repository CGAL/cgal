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

#ifndef CGAL_THREE_POINT_FAMILY_WEIGHTS_H
#define CGAL_THREE_POINT_FAMILY_WEIGHTS_H

#include <CGAL/Weights/internal/utils.h>

#include <CGAL/Point_2.h>
#include <CGAL/Point_3.h>

namespace CGAL {
namespace Weights {

/// \cond SKIP_IN_MANUAL
namespace three_point_family_ns {

template<typename FT>
FT weight(const FT d1, const FT d2, const FT d,
          const FT A1, const FT A2, const FT B,
          const FT p)
{
  FT w = FT(0);
  CGAL_precondition(!is_zero(A1) && !is_zero(A2));
  const FT prod = A1 * A2;
  if (!is_zero(prod))
  {
    const FT r1 = internal::power(d1, p);
    const FT r2 = internal::power(d2, p);
    const FT r  = internal::power(d , p);

    w = (r2 * A1 - r * B + r1 * A2) / prod;
  }
  return w;
}

} // namespace three_point_family_ns

template<typename GeomTraits>
typename GeomTraits::FT three_point_family_weight(const typename GeomTraits::Point_2& t,
                                                  const typename GeomTraits::Point_2& r,
                                                  const typename GeomTraits::Point_2& p,
                                                  const typename GeomTraits::Point_2& q,
                                                  const typename GeomTraits::FT a,
                                                  const GeomTraits& traits)
{
  using FT = typename GeomTraits::FT;

  const FT d1 = internal::distance_2(traits, q, t);
  const FT d2 = internal::distance_2(traits, q, r);
  const FT d3 = internal::distance_2(traits, q, p);

  const FT A1 = internal::area_2(traits, r, q, t);
  const FT A2 = internal::area_2(traits, p, q, r);
  const FT B  = internal::area_2(traits, p, q, t);

  return three_point_family_ns::weight(traits, d1, d2, d3, A1, A2, B, a);
}

template<typename GeomTraits>
typename GeomTraits::FT three_point_family_weight(const CGAL::Point_2<GeomTraits>& t,
                                                  const CGAL::Point_2<GeomTraits>& r,
                                                  const CGAL::Point_2<GeomTraits>& p,
                                                  const CGAL::Point_2<GeomTraits>& q,
                                                  const typename GeomTraits::FT a = typename GeomTraits::FT(1))
{
  const GeomTraits traits;
  return three_point_family_weight(t, r, p, q, a, traits);
}

namespace internal {

template<typename GeomTraits>
typename GeomTraits::FT three_point_family_weight(const typename GeomTraits::Point_3& t,
                                                  const typename GeomTraits::Point_3& r,
                                                  const typename GeomTraits::Point_3& p,
                                                  const typename GeomTraits::Point_3& q,
                                                  const typename GeomTraits::FT a,
                                                  const GeomTraits& traits)
{
  using Point_2 = typename GeomTraits::Point_2;

  Point_2 tf, rf, pf, qf;
  internal::flatten(traits,
                    t,  r,  p,  q,
                    tf, rf, pf, qf);
  return CGAL::Weights::three_point_family_weight(tf, rf, pf, qf, a, traits);
}

template<typename GeomTraits>
typename GeomTraits::FT three_point_family_weight(const CGAL::Point_3<GeomTraits>& t,
                                                  const CGAL::Point_3<GeomTraits>& r,
                                                  const CGAL::Point_3<GeomTraits>& p,
                                                  const CGAL::Point_3<GeomTraits>& q,
                                                  const typename GeomTraits::FT a = typename GeomTraits::FT(1))
{
  const GeomTraits traits;
  return three_point_family_weight(t, r, p, q, a, traits);
}

} // namespace internal

/// \endcond

} // namespace Weights
} // namespace CGAL

#endif // CGAL_THREE_POINT_FAMILY_WEIGHTS_H
