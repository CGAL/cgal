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
FT weight(const FT d0, const FT d2, const FT d,
          const FT A0, const FT A2, const FT B,
          const FT p)
{
  FT w = FT(0);
  CGAL_precondition(!is_zero(A0) && !is_zero(A2));
  const FT prod = A0 * A2;
  if (!is_zero(prod))
  {
    const FT r0 = internal::power(d0, p);
    const FT r = internal::power(d , p);
    const FT r2 = internal::power(d2, p);

    w = (r2 * A0 - r * B + r0 * A2) / prod;
  }
  return w;
}

} // namespace three_point_family_ns

/// \endcond

// 2D ==============================================================================================

/*!
  \ingroup PkgWeightsRefThreePointFamilyWeights
  \brief computes the three-point family weight in 2D at `q` using the points `p0`, `p1` and `p2`,
         and the power parameter `a`.
  \tparam GeomTraits a model of `AnalyticWeightTraits_2`
*/
template<typename GeomTraits>
typename GeomTraits::FT three_point_family_weight(const typename GeomTraits::Point_2& p0,
                                                  const typename GeomTraits::Point_2& p1,
                                                  const typename GeomTraits::Point_2& p2,
                                                  const typename GeomTraits::Point_2& q,
                                                  const typename GeomTraits::FT a,
                                                  const GeomTraits& traits)
{
  using FT = typename GeomTraits::FT;

  auto area_2 = traits.compute_area_2_object();

  const FT d0 = internal::distance_2(q, p0, traits);
  const FT d = internal::distance_2(q, p1, traits);
  const FT d2 = internal::distance_2(q, p2, traits);

  const FT A0 = area_2(p1, q, p0);
  const FT A2 = area_2(p2, q, p1);
  const FT B = area_2(p2, q, p0);

  return three_point_family_ns::weight(d0, d2, d, A0, A2, B, a);
}

/*!
  \ingroup PkgWeightsRefThreePointFamilyWeights
  \brief computes the three-point family weight in 2D at `q` using the points `p0`, `p1` and `p2`,
         and the power parameter `a`.
  \tparam Kernel a model of `Kernel`
*/
template<typename Kernel>
typename Kernel::FT three_point_family_weight(const CGAL::Point_2<Kernel>& p0,
                                              const CGAL::Point_2<Kernel>& p1,
                                              const CGAL::Point_2<Kernel>& p2,
                                              const CGAL::Point_2<Kernel>& q,
                                              const typename Kernel::FT a = {1})
{
  const Kernel traits;
  return three_point_family_weight(p0, p1, p2, q, a, traits);
}

// 3D ==============================================================================================

/// \cond SKIP_IN_MANUAL

template<typename GeomTraits>
typename GeomTraits::FT three_point_family_weight(const typename GeomTraits::Point_3& p0,
                                                  const typename GeomTraits::Point_3& p1,
                                                  const typename GeomTraits::Point_3& p2,
                                                  const typename GeomTraits::Point_3& q,
                                                  const typename GeomTraits::FT a,
                                                  const GeomTraits& traits)
{
  using Point_2 = typename GeomTraits::Point_2;

  Point_2 p0f, p1f, p2f, qf;
  internal::flatten(p0, p1, p2 , q,
                    p0f, p1f, p2f, qf,
                    traits);

  return CGAL::Weights::three_point_family_weight(p0f, p1f, p2f, qf, a, traits);
}

template<typename Kernel>
typename Kernel::FT three_point_family_weight(const CGAL::Point_3<Kernel>& p0,
                                              const CGAL::Point_3<Kernel>& p1,
                                              const CGAL::Point_3<Kernel>& p2,
                                              const CGAL::Point_3<Kernel>& q,
                                              const typename Kernel::FT a = {1})
{
  const Kernel traits;
  return three_point_family_weight(p0, p1, p2, q, a, traits);
}

/// \endcond

} // namespace Weights
} // namespace CGAL

#endif // CGAL_THREE_POINT_FAMILY_WEIGHTS_H
