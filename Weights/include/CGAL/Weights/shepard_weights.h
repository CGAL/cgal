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

template<typename FT>
FT weight(const FT d, const FT p)
{
  FT w = FT(0);
  CGAL_precondition(is_positive(d));
  if(is_positive(d))
    w = internal::power(d, -p);

  return w;
}

} // namespace shepard_ns

/// \endcond

// 2D ==============================================================================================

/*!
  \ingroup PkgWeightsRefShepardWeights
  \brief computes the Shepard weight in 2D using the points `p` and `q` and the power parameter `a`.
  \tparam GeomTraits a model of `AnalyticWeightTraits_2`
*/
template<typename GeomTraits>
typename GeomTraits::FT shepard_weight(const typename GeomTraits::Point_2&,
                                       const typename GeomTraits::Point_2& p,
                                       const typename GeomTraits::Point_2&,
                                       const typename GeomTraits::Point_2& q,
                                       const typename GeomTraits::FT a,
                                       const GeomTraits& traits)
{
  using FT = typename GeomTraits::FT;
  const FT d = internal::distance_2(p, q, traits);
  return shepard_ns::weight(d, a);
}

/*!
  \ingroup PkgWeightsRefShepardWeights
  \brief computes the Shepard weight in 2D using the points `p` and `q`, and the power parameter `a`
  \tparam Kernel a model of `Kernel`
*/
template<typename Kernel>
#ifdef DOXYGEN_RUNNING
typename Kernel::FT shepard_weight(const CGAL::Point_2<Kernel>&,
                                   const CGAL::Point_2<Kernel>& p^,
                                   const CGAL::Point_2<Kernel>&,
                                   const CGAL::Point_2<Kernel>& q,
                                   const typename Kernel::FT a = {1})
#else
typename Kernel::FT shepard_weight(const CGAL::Point_2<Kernel>& stub_l,
                                   const CGAL::Point_2<Kernel>& p,
                                   const CGAL::Point_2<Kernel>& stub_r,
                                   const CGAL::Point_2<Kernel>& q,
                                   const typename Kernel::FT a = {1})
#endif
{
  const Kernel traits;
  return shepard_weight(stub_l, p, stub_r, q, a, traits);
}

/*!
  \ingroup PkgWeightsRefShepardWeights
  \brief computes the Shepard weight in 2D using the points `p` and `q` and the power parameter `a`.
  \tparam GeomTraits a model of `AnalyticWeightTraits_2`
*/
template<typename GeomTraits>
typename GeomTraits::FT shepard_weight(const typename GeomTraits::Point_2& p,
                                       const typename GeomTraits::Point_2& q,
                                       const typename GeomTraits::FT a,
                                       const GeomTraits& traits)
{
  typename GeomTraits::Point_2 stub;
  return shepard_weight(stub, p, stub, q, a, traits);
}

/*!
  \ingroup PkgWeightsRefShepardWeights
  \brief computes the Shepard weight in 2D using the points `p` and `q`, and the power parameter `a`.
  \tparam Kernel a model of `Kernel`
*/
template<typename Kernel>
typename Kernel::FT shepard_weight(const CGAL::Point_2<Kernel>& p,
                                   const CGAL::Point_2<Kernel>& q,
                                   const typename Kernel::FT a = {1})
{
  CGAL::Point_2<Kernel> stub;
  return shepard_weight(stub, p, stub, q, a);
}

// 3D ==============================================================================================

/*!
  \ingroup PkgWeightsRefShepardWeights
  \brief computes the Shepard weight in 3D using the points `p` and `q` and the power parameter `a`.
  \tparam GeomTraits a model of `AnalyticWeightTraits_3`
*/
template<typename GeomTraits>
typename GeomTraits::FT shepard_weight(const typename GeomTraits::Point_3&,
                                       const typename GeomTraits::Point_3& p,
                                       const typename GeomTraits::Point_3&,
                                       const typename GeomTraits::Point_3& q,
                                       const typename GeomTraits::FT a,
                                       const GeomTraits& traits)
{
  using FT = typename GeomTraits::FT;
  const FT d = internal::distance_3(p, q, traits);
  return shepard_ns::weight(d, a);
}

/*!
  \ingroup PkgWeightsRefShepardWeights
  \brief computes the Shepard weight in 3D using the points `p` and `q`, and the power parameter `a`.
  \tparam Kernel a model of `Kernel`
*/
template<typename Kernel>
#ifdef DOXYGEN_RUNNING
typename Kernel::FT shepard_weight(const CGAL::Point_3<Kernel>& p,
                                   const CGAL::Point_3<Kernel>&,
                                   const CGAL::Point_3<Kernel>& q,
                                   const CGAL::Point_3<Kernel>&,
                                   const typename Kernel::FT a = {1})
#else
typename Kernel::FT shepard_weight(const CGAL::Point_3<Kernel>& stub_l,
                                   const CGAL::Point_3<Kernel>& p,
                                   const CGAL::Point_3<Kernel>& stub_r,
                                   const CGAL::Point_3<Kernel>& q,
                                   const typename Kernel::FT a = {1})
#endif
{
  const Kernel traits;
  return shepard_weight(stub_l, p, stub_r, q, a, traits);
}

/*!
  \ingroup PkgWeightsRefShepardWeights
  \brief computes the Shepard weight in 3D using the points `p` and `q` and the power parameter `a`.
  \tparam GeomTraits a model of `AnalyticWeightTraits_3`
*/
template<typename GeomTraits>
typename GeomTraits::FT shepard_weight(const typename GeomTraits::Point_3& p,
                                       const typename GeomTraits::Point_3& q,
                                       const typename GeomTraits::FT a,
                                       const GeomTraits& traits)
{
  typename GeomTraits::Point_3 stub;
  return shepard_weight(stub, p, stub, q, a, traits);
}

/*!
  \ingroup PkgWeightsRefShepardWeights
  \brief computes the Shepard weight in 3D using the points `p` and `q`, and the power parameter `a`.
  \tparam Kernel a model of `Kernel`
*/
template<typename Kernel>
typename Kernel::FT shepard_weight(const CGAL::Point_3<Kernel>& p,
                                   const CGAL::Point_3<Kernel>& q,
                                   const typename Kernel::FT a = {1})
{
  CGAL::Point_3<Kernel> stub;
  return shepard_weight(stub, p, stub, q, a);
}

} // namespace Weights
} // namespace CGAL

#endif // CGAL_SHEPARD_WEIGHTS_H
