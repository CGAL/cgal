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

// Internal includes.
#include <CGAL/Weights/internal/utils.h>

namespace CGAL {
namespace Weights {

  /// \cond SKIP_IN_MANUAL
  namespace shepard_ns {

    template<typename GeomTraits>
    typename GeomTraits::FT weight(
      const GeomTraits& traits,
      const typename GeomTraits::FT d,
      const typename GeomTraits::FT p) {

      using FT = typename GeomTraits::FT;
      FT w = FT(0);
      CGAL_precondition(d != FT(0));
      if (d != FT(0)) {
        FT denom = d;
        if (p != FT(1)) {
          denom = internal::power(traits, d, p);
        }
        w = FT(1) / denom;
      }
      return w;
    }
  }
  /// \endcond

  #if defined(DOXYGEN_RUNNING)

  /*!
    \ingroup PkgWeightsRefShepardWeights

    \brief computes the Shepard weight in 2D using the points `p` and `q` and the power parameter `a`,
    given a traits class `traits` with geometric objects, predicates, and constructions.
  */
  template<typename GeomTraits>
  typename GeomTraits::FT shepard_weight(
    const typename GeomTraits::Point_2&,
    const typename GeomTraits::Point_2& p,
    const typename GeomTraits::Point_2&,
    const typename GeomTraits::Point_2& q,
    const typename GeomTraits::FT a,
    const GeomTraits& traits) { }

  /*!
    \ingroup PkgWeightsRefShepardWeights

    \brief computes the Shepard weight in 3D using the points `p` and `q` and the power parameter `a`,
    given a traits class `traits` with geometric objects, predicates, and constructions.
  */
  template<typename GeomTraits>
  typename GeomTraits::FT shepard_weight(
    const typename GeomTraits::Point_3&,
    const typename GeomTraits::Point_3& p,
    const typename GeomTraits::Point_3&,
    const typename GeomTraits::Point_3& q,
    const typename GeomTraits::FT a,
    const GeomTraits& traits) { }

  /*!
    \ingroup PkgWeightsRefShepardWeights

    \brief computes the Shepard weight in 2D using the points `p` and `q`,
    which are parameterized by a `Kernel` K, and the power parameter `a` which
    can be omitted.
  */
  template<typename K>
  typename K::FT shepard_weight(
    const CGAL::Point_2<K>&,
    const CGAL::Point_2<K>& p,
    const CGAL::Point_2<K>&,
    const CGAL::Point_2<K>& q,
    const typename K::FT a = typename K::FT(1)) { }

  /*!
    \ingroup PkgWeightsRefShepardWeights

    \brief computes the Shepard weight in 3D using the points `p` and `q`,
    which are parameterized by a `Kernel` K, and the power parameter `a` which
    can be omitted.
  */
  template<typename K>
  typename K::FT shepard_weight(
    const CGAL::Point_3<K>&,
    const CGAL::Point_3<K>& p,
    const CGAL::Point_3<K>&,
    const CGAL::Point_3<K>& q,
    const typename K::FT a = typename K::FT(1)) { }

  /*!
    \ingroup PkgWeightsRefShepardWeights

    \brief computes the Shepard weight in 2D using the points `p` and `q` and the power parameter `a`,
    given a traits class `traits` with geometric objects, predicates, and constructions.
  */
  template<typename GeomTraits>
  typename GeomTraits::FT shepard_weight(
    const typename GeomTraits::Point_2& p,
    const typename GeomTraits::Point_2& q,
    const typename GeomTraits::FT a,
    const GeomTraits& traits) { }

  /*!
    \ingroup PkgWeightsRefShepardWeights

    \brief computes the Shepard weight in 3D using the points `p` and `q` and the power parameter `a`,
    given a traits class `traits` with geometric objects, predicates, and constructions.
  */
  template<typename GeomTraits>
  typename GeomTraits::FT shepard_weight(
    const typename GeomTraits::Point_3& p,
    const typename GeomTraits::Point_3& q,
    const typename GeomTraits::FT a,
    const GeomTraits& traits) { }

  /*!
    \ingroup PkgWeightsRefShepardWeights

    \brief computes the Shepard weight in 2D using the points `p` and `q`,
    which are parameterized by a `Kernel` K, and the power parameter `a` which
    can be omitted.
  */
  template<typename K>
  typename K::FT shepard_weight(
    const CGAL::Point_2<K>& p,
    const CGAL::Point_2<K>& q,
    const typename K::FT a = typename K::FT(1)) { }

  /*!
    \ingroup PkgWeightsRefShepardWeights

    \brief computes the Shepard weight in 3D using the points `p` and `q`,
    which are parameterized by a `Kernel` K, and the power parameter `a` which
    can be omitted.
  */
  template<typename K>
  typename K::FT shepard_weight(
    const CGAL::Point_3<K>& p,
    const CGAL::Point_3<K>& q,
    const typename K::FT a = typename K::FT(1)) { }

  #endif // DOXYGEN_RUNNING

  /// \cond SKIP_IN_MANUAL
  template<typename GeomTraits>
  typename GeomTraits::FT shepard_weight(
    const typename GeomTraits::Point_2&,
    const typename GeomTraits::Point_2& r,
    const typename GeomTraits::Point_2&,
    const typename GeomTraits::Point_2& q,
    const typename GeomTraits::FT a,
    const GeomTraits& traits) {

    using FT = typename GeomTraits::FT;
    const FT d = internal::distance_2(traits, q, r);
    return shepard_ns::weight(traits, d, a);
  }

  template<typename GeomTraits>
  typename GeomTraits::FT shepard_weight(
    const CGAL::Point_2<GeomTraits>& t,
    const CGAL::Point_2<GeomTraits>& r,
    const CGAL::Point_2<GeomTraits>& p,
    const CGAL::Point_2<GeomTraits>& q,
    const typename GeomTraits::FT a =
    typename GeomTraits::FT(1)) {

    const GeomTraits traits;
    return shepard_weight(t, r, p, q, a, traits);
  }

  template<typename GeomTraits>
  typename GeomTraits::FT shepard_weight(
    const typename GeomTraits::Point_2& p,
    const typename GeomTraits::Point_2& q,
    const typename GeomTraits::FT a,
    const GeomTraits& traits) {

    typename GeomTraits::Point_2 stub;
    return shepard_weight(stub, p, stub, q, a, traits);
  }

  template<typename GeomTraits>
  typename GeomTraits::FT shepard_weight(
    const CGAL::Point_2<GeomTraits>& p,
    const CGAL::Point_2<GeomTraits>& q,
    const typename GeomTraits::FT a =
    typename GeomTraits::FT(1)) {

    CGAL::Point_2<GeomTraits> stub;
    return shepard_weight(stub, p, stub, q, a);
  }

  template<typename GeomTraits>
  typename GeomTraits::FT shepard_weight(
    const typename GeomTraits::Point_3&,
    const typename GeomTraits::Point_3& r,
    const typename GeomTraits::Point_3&,
    const typename GeomTraits::Point_3& q,
    const typename GeomTraits::FT a,
    const GeomTraits& traits) {

    using FT = typename GeomTraits::FT;
    const FT d = internal::distance_3(traits, q, r);
    return shepard_ns::weight(traits, d, a);
  }

  template<typename GeomTraits>
  typename GeomTraits::FT shepard_weight(
    const CGAL::Point_3<GeomTraits>& t,
    const CGAL::Point_3<GeomTraits>& r,
    const CGAL::Point_3<GeomTraits>& p,
    const CGAL::Point_3<GeomTraits>& q,
    const typename GeomTraits::FT a =
    typename GeomTraits::FT(1)) {

    const GeomTraits traits;
    return shepard_weight(t, r, p, q, a, traits);
  }

  template<typename GeomTraits>
  typename GeomTraits::FT shepard_weight(
    const typename GeomTraits::Point_3& p,
    const typename GeomTraits::Point_3& q,
    const typename GeomTraits::FT a,
    const GeomTraits& traits) {

    typename GeomTraits::Point_3 stub;
    return shepard_weight(stub, p, stub, q, a, traits);
  }

  template<typename GeomTraits>
  typename GeomTraits::FT shepard_weight(
    const CGAL::Point_3<GeomTraits>& p,
    const CGAL::Point_3<GeomTraits>& q,
    const typename GeomTraits::FT a =
    typename GeomTraits::FT(1)) {

    CGAL::Point_3<GeomTraits> stub;
    return shepard_weight(stub, p, stub, q, a);
  }
  /// \endcond

} // namespace Weights
} // namespace CGAL

#endif // CGAL_SHEPARD_WEIGHTS_H
