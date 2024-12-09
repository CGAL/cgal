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
  CGAL_precondition(!is_zero(d));
  if (!is_zero(d))
    w = FT(1) / d;

  return w;
}

} // namespace inverse_distance_ns

/// \endcond

/*!
  \ingroup PkgWeightsRefInverseDistanceWeights
  \brief computes the inverse distance weight in 2D using the points `p` and `q`.
  \tparam GeomTraits a model of `AnalyticWeightTraits_2`
*/
template<typename GeomTraits>
typename GeomTraits::FT inverse_distance_weight(const typename GeomTraits::Point_2&,
                                                const typename GeomTraits::Point_2& p,
                                                const typename GeomTraits::Point_2&,
                                                const typename GeomTraits::Point_2& q,
                                                const GeomTraits& traits)
{
  using FT = typename GeomTraits::FT;

  const FT d = internal::distance_2(p, q, traits);
  return inverse_distance_ns::weight(d);
}

/*!
  \ingroup PkgWeightsRefInverseDistanceWeights
  \brief computes the inverse distance weight in 2D using the points `p` and `q`.
  \tparam Kernel a model of `Kernel`
*/
template<typename Kernel>
#ifdef DOXYGEN_RUNNING
typename Kernel::FT inverse_distance_weight(const CGAL::Point_2<Kernel>&,
                                            const CGAL::Point_2<Kernel>& p,
                                            const CGAL::Point_2<Kernel>&,
                                            const CGAL::Point_2<Kernel>& q)
#else
typename Kernel::FT inverse_distance_weight(const CGAL::Point_2<Kernel>& stub_l,
                                            const CGAL::Point_2<Kernel>& p,
                                            const CGAL::Point_2<Kernel>& stub_r,
                                            const CGAL::Point_2<Kernel>& q)
#endif
{
  const Kernel traits;
  return inverse_distance_weight(stub_l, p, stub_r, q, traits);
}

/*!
  \ingroup PkgWeightsRefInverseDistanceWeights
  \brief computes the inverse distance weight in 2D using the points `p` and `q`.
  \tparam GeomTraits a model of `AnalyticWeightTraits_2`
*/
template<typename GeomTraits>
typename GeomTraits::FT inverse_distance_weight(const typename GeomTraits::Point_2& p,
                                                const typename GeomTraits::Point_2& q,
                                                const GeomTraits& traits)
{
  typename GeomTraits::Point_2 stub;
  return inverse_distance_weight(stub, p, stub, q, traits);
}

/*!
  \ingroup PkgWeightsRefInverseDistanceWeights
  \brief computes the inverse distance weight in 2D using the points `p` and `q`.
  \tparam Kernel a model of `Kernel`
*/
template<typename Kernel>
typename Kernel::FT inverse_distance_weight(const CGAL::Point_2<Kernel>& p,
                                            const CGAL::Point_2<Kernel>& q)
{
  CGAL::Point_2<Kernel> stub;
  return inverse_distance_weight(stub, p, stub, q);
}

// 3D ==============================================================================================

/*!
  \ingroup PkgWeightsRefInverseDistanceWeights
  \brief computes the inverse distance weight in 3D using the points `p` and `q`.
  \tparam GeomTraits a model of `AnalyticWeightTraits_3`
*/
template<typename GeomTraits>
typename GeomTraits::FT inverse_distance_weight(const typename GeomTraits::Point_3&,
                                                const typename GeomTraits::Point_3& p,
                                                const typename GeomTraits::Point_3&,
                                                const typename GeomTraits::Point_3& q,
                                                const GeomTraits& traits)
{
  using FT = typename GeomTraits::FT;

  const FT d = internal::distance_3(p, q, traits);
  return inverse_distance_ns::weight(d);
}

/*!
  \ingroup PkgWeightsRefInverseDistanceWeights
  \brief computes the inverse distance weight in 3D using the points `p` and `q`.
  \tparam Kernel a model of `Kernel`
*/
template<typename Kernel>
#ifdef DOXYGEN_RUNNING
typename Kernel::FT inverse_distance_weight(const CGAL::Point_3<Kernel>&,
                                            const CGAL::Point_3<Kernel>& p,
                                            const CGAL::Point_3<Kernel>&,
                                            const CGAL::Point_3<Kernel>& q)
#else
typename Kernel::FT inverse_distance_weight(const CGAL::Point_3<Kernel>& stub_l,
                                            const CGAL::Point_3<Kernel>& p,
                                            const CGAL::Point_3<Kernel>& stub_r,
                                            const CGAL::Point_3<Kernel>& q)
#endif
{
  const Kernel traits;
  return inverse_distance_weight(stub_l, p, stub_r, q, traits);
}

/*!
  \ingroup PkgWeightsRefInverseDistanceWeights
  \brief computes the inverse distance weight in 3D using the points `p` and `q`.
  \tparam GeomTraits a model of `AnalyticWeightTraits_3`
*/
template<typename GeomTraits>
typename GeomTraits::FT inverse_distance_weight(const typename GeomTraits::Point_3& p,
                                                const typename GeomTraits::Point_3& q,
                                                const GeomTraits& traits)
{
  typename GeomTraits::Point_3 stub;
  return inverse_distance_weight(stub, p, stub, q, traits);
}

/*!
  \ingroup PkgWeightsRefInverseDistanceWeights
  \brief computes the inverse distance weight in 3D using the points `p` and `q`.
  \tparam Kernel a model of `Kernel`
*/
template<typename Kernel>
typename Kernel::FT inverse_distance_weight(const CGAL::Point_3<Kernel>& p,
                                            const CGAL::Point_3<Kernel>& q)
{
  CGAL::Point_3<Kernel> stub;
  return inverse_distance_weight(stub, p, stub, q);
}

} // namespace Weights
} // namespace CGAL

#endif // CGAL_INVERSE_DISTANCE_WEIGHTS_H
