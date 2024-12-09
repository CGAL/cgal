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

#ifndef CGAL_AUTHALIC_WEIGHTS_H
#define CGAL_AUTHALIC_WEIGHTS_H

#include <CGAL/Weights/utils.h>

#include <CGAL/Point_2.h>
#include <CGAL/Point_3.h>

namespace CGAL {
namespace Weights {

/// \cond SKIP_IN_MANUAL

namespace authalic_ns {

template<typename FT>
FT half_weight(const FT cot, const FT r2)
{
  FT w = FT(0);
  CGAL_precondition(!is_zero(r2));
  if (!is_zero(r2))
    w = FT(2) * cot / r2;

  return w;
}

template<typename FT>
FT weight(const FT cot_gamma, const FT cot_beta, const FT r2)
{
  FT w = FT(0);
  CGAL_precondition(!is_zero(r2));
  if (!is_zero(r2))
    w = FT(2) * (cot_gamma + cot_beta) / r2;

  return w;
}

} // namespace authalic_ns

/// \endcond

/*!
  \ingroup PkgWeightsRefAuthalicWeights

  \brief computes the half value of the authalic weight.

  This function computes the half of the authalic weight using the precomputed
  cotangent and squared distance values. The returned value is
  \f$\frac{2\textbf{cot}}{\textbf{sq_d}}\f$.

  \tparam FT a model of `FieldNumberType`

  \param cot the cotangent value
  \param sq_d the squared distance value

  \pre sq_d != 0

  \sa `authalic_weight()`
*/
template<typename FT>
FT half_authalic_weight(const FT cot, const FT sq_d)
{
  return authalic_ns::half_weight(cot, sq_d);
}

// 2D ==============================================================================================

/*!
  \ingroup PkgWeightsRefAuthalicWeights
  \brief computes the authalic weight in 2D at `q` using the points `p0`, `p1`, and `p2`.
  \tparam GeomTraits a model of `AnalyticWeightTraits_2`
*/
template<typename GeomTraits>
typename GeomTraits::FT authalic_weight(const typename GeomTraits::Point_2& p0,
                                        const typename GeomTraits::Point_2& p1,
                                        const typename GeomTraits::Point_2& p2,
                                        const typename GeomTraits::Point_2& q,
                                        const GeomTraits& traits)
{
  using FT = typename GeomTraits::FT;

  const FT cot_gamma = cotangent_2(p0, p1, q, traits);
  const FT cot_beta = cotangent_2(q, p1, p2, traits);

  auto squared_distance_2 = traits.compute_squared_distance_2_object();
  const FT sq_d = squared_distance_2(q, p1);

  return authalic_ns::weight(cot_gamma, cot_beta, sq_d);
}

/*!
  \ingroup PkgWeightsRefAuthalicWeights
  \brief computes the authalic weight in 2D at `q` using the points `p0`, `p1`, and `p2`.
  \tparam Kernel a model of `Kernel`
*/
template<typename Kernel>
typename Kernel::FT authalic_weight(const CGAL::Point_2<Kernel>& p0,
                                    const CGAL::Point_2<Kernel>& p1,
                                    const CGAL::Point_2<Kernel>& p2,
                                    const CGAL::Point_2<Kernel>& q)
{
  const Kernel traits;
  return authalic_weight(p0, p1, p2, q, traits);
}

// 3D ==============================================================================================

/*!
  \ingroup PkgWeightsRefAuthalicWeights
  \brief computes the authalic weight in 3D at `q` using the points `p0`, `p1`, and `p2`.
  \tparam GeomTraits a model of `AnalyticWeightTraits_3`
*/
template<typename GeomTraits>
typename GeomTraits::FT authalic_weight(const typename GeomTraits::Point_3& p0,
                                        const typename GeomTraits::Point_3& p1,
                                        const typename GeomTraits::Point_3& p2,
                                        const typename GeomTraits::Point_3& q,
                                        const GeomTraits& traits)
{
  using FT = typename GeomTraits::FT;

  const FT cot_gamma = cotangent_3(p0, p1, q, traits);
  const FT cot_beta  = cotangent_3(q, p1, p2, traits);

  auto squared_distance_3 = traits.compute_squared_distance_3_object();
  const FT sq_d = squared_distance_3(q, p1);

  return authalic_ns::weight(cot_gamma, cot_beta, sq_d);
}

/*!
  \ingroup PkgWeightsRefAuthalicWeights
  \brief computes the authalic weight in 3D at `q` using the points `p0`, `p1`, and `p2`.
  \tparam Kernel a model of `Kernel`
*/
template<typename Kernel>
typename Kernel::FT authalic_weight(const CGAL::Point_3<Kernel>& p0,
                                    const CGAL::Point_3<Kernel>& p1,
                                    const CGAL::Point_3<Kernel>& p2,
                                    const CGAL::Point_3<Kernel>& q)
{
  const Kernel traits;
  return authalic_weight(p0, p1, p2, q, traits);
}

} // namespace Weights
} // namespace CGAL

#endif // CGAL_AUTHALIC_WEIGHTS_H
