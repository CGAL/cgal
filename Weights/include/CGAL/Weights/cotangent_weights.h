// Copyright (c) 2020 GeometryFactory SARL (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Dmitry Anisimov
//

#ifndef CGAL_COTANGENT_WEIGHTS_H
#define CGAL_COTANGENT_WEIGHTS_H

// #include <CGAL/license/Weights.h>

// Internal includes.
#include <CGAL/Weights/internal/utils.h>

namespace CGAL {
namespace Weights {

  /// \cond SKIP_IN_MANUAL
  namespace cotangent_ns {

    template<typename FT>
    FT half_weight(const FT cot) {
      return FT(2) * cot;
    }

    template<typename FT>
    FT weight(const FT cot_beta, const FT cot_gamma) {
      return FT(2) * (cot_beta + cot_gamma);
    }
  }
  /// \endcond

  /*!
    \ingroup PkgWeightsRefCotangentWeights

    \brief computes the half value of the cotangent weight.

    This function constructs the half of the cotangent weight using the precomputed
    cotangent value. The returned value is
    \f$2\textbf{cot}\f$.

    \tparam FT
    a model of `FieldNumberType`

    \param cot
    the cotangent value

    \sa `cotangent_weight()`
  */
  template<typename FT>
  FT half_cotangent_weight(const FT cot) {
    return cotangent_ns::half_weight(cot);
  }

  #if defined(DOXYGEN_RUNNING)

  /*!
    \ingroup PkgWeightsRefCotangentWeights

    \brief computes the cotangent weight in 2D at `q` using the points `p0`, `p1`,
    and `p2`, given a traits class `traits` with geometric objects, predicates, and constructions.
  */
  template<typename GeomTraits>
  typename GeomTraits::FT cotangent_weight(
    const typename GeomTraits::Point_2& p0,
    const typename GeomTraits::Point_2& p1,
    const typename GeomTraits::Point_2& p2,
    const typename GeomTraits::Point_2& q,
    const GeomTraits& traits) { }

  /*!
    \ingroup PkgWeightsRefCotangentWeights

    \brief computes the cotangent weight in 3D at `q` using the points `p0`, `p1`,
    and `p2`, given a traits class `traits` with geometric objects, predicates, and constructions.
  */
  template<typename GeomTraits>
  typename GeomTraits::FT cotangent_weight(
    const typename GeomTraits::Point_3& p0,
    const typename GeomTraits::Point_3& p1,
    const typename GeomTraits::Point_3& p2,
    const typename GeomTraits::Point_3& q,
    const GeomTraits& traits) { }

  /*!
    \ingroup PkgWeightsRefCotangentWeights

    \brief computes the cotangent weight in 2D at `q` using the points `p0`, `p1`,
    and `p2` which are parameterized by a `Kernel` K.
  */
  template<typename K>
  typename K::FT cotangent_weight(
    const CGAL::Point_2<K>& p0,
    const CGAL::Point_2<K>& p1,
    const CGAL::Point_2<K>& p2,
    const CGAL::Point_2<K>& q) { }

  /*!
    \ingroup PkgWeightsRefCotangentWeights

    \brief computes the cotangent weight in 3D at `q` using the points `p0`, `p1`,
    and `p2` which are parameterized by a `Kernel` K.
  */
  template<typename K>
  typename K::FT cotangent_weight(
    const CGAL::Point_3<K>& p0,
    const CGAL::Point_3<K>& p1,
    const CGAL::Point_3<K>& p2,
    const CGAL::Point_3<K>& q) { }

  #endif // DOXYGEN_RUNNING

  /// \cond SKIP_IN_MANUAL
  template<typename GeomTraits>
  typename GeomTraits::FT cotangent_weight(
    const typename GeomTraits::Point_2& t,
    const typename GeomTraits::Point_2& r,
    const typename GeomTraits::Point_2& p,
    const typename GeomTraits::Point_2& q,
    const GeomTraits& traits) {

    using FT = typename GeomTraits::FT;
    const FT cot_beta  = internal::cotangent_2(traits, q, t, r);
    const FT cot_gamma = internal::cotangent_2(traits, r, p, q);
    return cotangent_ns::weight(cot_beta, cot_gamma);
  }

  template<typename GeomTraits>
  typename GeomTraits::FT cotangent_weight(
    const CGAL::Point_2<GeomTraits>& t,
    const CGAL::Point_2<GeomTraits>& r,
    const CGAL::Point_2<GeomTraits>& p,
    const CGAL::Point_2<GeomTraits>& q) {

    const GeomTraits traits;
    return cotangent_weight(t, r, p, q, traits);
  }

  template<typename GeomTraits>
  typename GeomTraits::FT cotangent_weight(
    const typename GeomTraits::Point_3& t,
    const typename GeomTraits::Point_3& r,
    const typename GeomTraits::Point_3& p,
    const typename GeomTraits::Point_3& q,
    const GeomTraits& traits) {

    using FT = typename GeomTraits::FT;
    const FT cot_beta  = internal::cotangent_3(traits, q, t, r);
    const FT cot_gamma = internal::cotangent_3(traits, r, p, q);
    return cotangent_ns::weight(cot_beta, cot_gamma);
  }

  template<typename GeomTraits>
  typename GeomTraits::FT cotangent_weight(
    const CGAL::Point_3<GeomTraits>& t,
    const CGAL::Point_3<GeomTraits>& r,
    const CGAL::Point_3<GeomTraits>& p,
    const CGAL::Point_3<GeomTraits>& q) {

    const GeomTraits traits;
    return cotangent_weight(t, r, p, q, traits);
  }
  /// \endcond

} // namespace Weights
} // namespace CGAL

#endif // CGAL_COTANGENT_WEIGHTS_H
