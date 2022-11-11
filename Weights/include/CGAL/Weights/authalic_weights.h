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

// Internal includes.
#include <CGAL/Weights/internal/utils.h>

namespace CGAL {
namespace Weights {

  /// \cond SKIP_IN_MANUAL
  namespace authalic_ns {

    template<typename FT>
    FT half_weight(const FT cot, const FT r2) {

      FT w = FT(0);
      CGAL_precondition(r2 != FT(0));
      if (r2 != FT(0)) {
        const FT inv = FT(2) / r2;
        w = cot * inv;
      }
      return w;
    }

    template<typename FT>
    FT weight(const FT cot_gamma, const FT cot_beta, const FT r2) {

      FT w = FT(0);
      CGAL_precondition(r2 != FT(0));
      if (r2 != FT(0)) {
        const FT inv = FT(2) / r2;
        w = (cot_gamma + cot_beta) * inv;
      }
      return w;
    }
  }
  /// \endcond

  /*!
    \ingroup PkgWeightsRefAuthalicWeights

    \brief computes the half value of the authalic weight.

    This function constructs the half of the authalic weight using the precomputed
    cotangent and squared distance values. The returned value is
    \f$\frac{2\textbf{cot}}{\textbf{d2}}\f$.

    \tparam FT
    a model of `FieldNumberType`

    \param cot
    the cotangent value

    \param d2
    the squared distance value

    \pre d2 != 0

    \sa `authalic_weight()`
  */
  template<typename FT>
  FT half_authalic_weight(const FT cot, const FT d2) {
    return authalic_ns::half_weight(cot, d2);
  }

  #if defined(DOXYGEN_RUNNING)

  /*!
    \ingroup PkgWeightsRefAuthalicWeights

    \brief computes the authalic weight in 2D at `q` using the points `p0`, `p1`,
    and `p2`, given a traits class `traits` with geometric objects, predicates, and constructions.
  */
  template<typename GeomTraits>
  typename GeomTraits::FT authalic_weight(
    const typename GeomTraits::Point_2& p0,
    const typename GeomTraits::Point_2& p1,
    const typename GeomTraits::Point_2& p2,
    const typename GeomTraits::Point_2& q,
    const GeomTraits& traits) { }

  /*!
    \ingroup PkgWeightsRefAuthalicWeights

    \brief computes the authalic weight in 3D at `q` using the points `p0`, `p1`,
    and `p2`, given a traits class `traits` with geometric objects, predicates, and constructions.
  */
  template<typename GeomTraits>
  typename GeomTraits::FT authalic_weight(
    const typename GeomTraits::Point_3& p0,
    const typename GeomTraits::Point_3& p1,
    const typename GeomTraits::Point_3& p2,
    const typename GeomTraits::Point_3& q,
    const GeomTraits& traits) { }

  /*!
    \ingroup PkgWeightsRefAuthalicWeights

    \brief computes the authalic weight in 2D at `q` using the points `p0`, `p1`,
    and `p2` which are parameterized by a `Kernel` K.
  */
  template<typename K>
  typename K::FT authalic_weight(
    const CGAL::Point_2<K>& p0,
    const CGAL::Point_2<K>& p1,
    const CGAL::Point_2<K>& p2,
    const CGAL::Point_2<K>& q) { }

  /*!
    \ingroup PkgWeightsRefAuthalicWeights

    \brief computes the authalic weight in 3D at `q` using the points `p0`, `p1`,
    and `p2` which are parameterized by a `Kernel` K.
  */
  template<typename K>
  typename K::FT authalic_weight(
    const CGAL::Point_3<K>& p0,
    const CGAL::Point_3<K>& p1,
    const CGAL::Point_3<K>& p2,
    const CGAL::Point_3<K>& q) { }

  #endif // DOXYGEN_RUNNING

  /// \cond SKIP_IN_MANUAL
  // Overloads!
  template<typename GeomTraits>
  typename GeomTraits::FT authalic_weight(
    const typename GeomTraits::Point_2& t,
    const typename GeomTraits::Point_2& r,
    const typename GeomTraits::Point_2& p,
    const typename GeomTraits::Point_2& q,
    const GeomTraits& traits) {

    using FT = typename GeomTraits::FT;
    const FT cot_gamma = internal::cotangent_2(traits, t, r, q);
    const FT cot_beta  = internal::cotangent_2(traits, q, r, p);

    const auto squared_distance_2 =
      traits.compute_squared_distance_2_object();
    const FT d2 = squared_distance_2(q, r);
    return authalic_ns::weight(cot_gamma, cot_beta, d2);
  }

  template<typename GeomTraits>
  typename GeomTraits::FT authalic_weight(
    const CGAL::Point_2<GeomTraits>& t,
    const CGAL::Point_2<GeomTraits>& r,
    const CGAL::Point_2<GeomTraits>& p,
    const CGAL::Point_2<GeomTraits>& q) {

    const GeomTraits traits;
    return authalic_weight(t, r, p, q, traits);
  }

  template<typename GeomTraits>
  typename GeomTraits::FT authalic_weight(
    const typename GeomTraits::Point_3& t,
    const typename GeomTraits::Point_3& r,
    const typename GeomTraits::Point_3& p,
    const typename GeomTraits::Point_3& q,
    const GeomTraits& traits) {

    using FT = typename GeomTraits::FT;
    const FT cot_gamma = internal::cotangent_3(traits, t, r, q);
    const FT cot_beta  = internal::cotangent_3(traits, q, r, p);

    const auto squared_distance_3 =
      traits.compute_squared_distance_3_object();
    const FT d2 = squared_distance_3(q, r);
    return authalic_ns::weight(cot_gamma, cot_beta, d2);
  }

  template<typename GeomTraits>
  typename GeomTraits::FT authalic_weight(
    const CGAL::Point_3<GeomTraits>& t,
    const CGAL::Point_3<GeomTraits>& r,
    const CGAL::Point_3<GeomTraits>& p,
    const CGAL::Point_3<GeomTraits>& q) {

    const GeomTraits traits;
    return authalic_weight(t, r, p, q, traits);
  }
  /// \endcond

} // namespace Weights
} // namespace CGAL

#endif // CGAL_AUTHALIC_WEIGHTS_H
