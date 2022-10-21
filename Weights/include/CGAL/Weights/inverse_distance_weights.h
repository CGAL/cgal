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

// Internal includes.
#include <CGAL/Weights/internal/utils.h>

namespace CGAL {
namespace Weights {

  /// \cond SKIP_IN_MANUAL
  namespace inverse_distance_ns {

    template<typename FT>
    FT weight(const FT d) {

      FT w = FT(0);
      CGAL_precondition(d != FT(0));
      if (d != FT(0)) {
        w = FT(1) / d;
      }
      return w;
    }
  }
  /// \endcond

  #if defined(DOXYGEN_RUNNING)

  /*!
    \ingroup PkgWeightsRefInverseDistanceWeights

    \brief computes the inverse distance weight in 2D using the points `p` and `q`,
    given a traits class `traits` with geometric objects, predicates, and constructions.
  */
  template<typename GeomTraits>
  typename GeomTraits::FT inverse_distance_weight(
    const typename GeomTraits::Point_2&,
    const typename GeomTraits::Point_2& p,
    const typename GeomTraits::Point_2&,
    const typename GeomTraits::Point_2& q,
    const GeomTraits& traits) { }

  /*!
    \ingroup PkgWeightsRefInverseDistanceWeights

    \brief computes the inverse distance weight in 3D using the points `p` and `q`,
    given a traits class `traits` with geometric objects, predicates, and constructions.
  */
  template<typename GeomTraits>
  typename GeomTraits::FT inverse_distance_weight(
    const typename GeomTraits::Point_3&,
    const typename GeomTraits::Point_3& p,
    const typename GeomTraits::Point_3&,
    const typename GeomTraits::Point_3& q,
    const GeomTraits& traits) { }

  /*!
    \ingroup PkgWeightsRefInverseDistanceWeights

    \brief computes the inverse distance weight in 2D using the points `p` and `q`,
    which are parameterized by a `Kernel` K.
  */
  template<typename K>
  typename K::FT inverse_distance_weight(
    const CGAL::Point_2<K>&,
    const CGAL::Point_2<K>& p,
    const CGAL::Point_2<K>&,
    const CGAL::Point_2<K>& q) { }

  /*!
    \ingroup PkgWeightsRefInverseDistanceWeights

    \brief computes the inverse distance weight in 3D using the points `p` and `q`,
    which are parameterized by a `Kernel` K.
  */
  template<typename K>
  typename K::FT inverse_distance_weight(
    const CGAL::Point_3<K>&,
    const CGAL::Point_3<K>& p,
    const CGAL::Point_3<K>&,
    const CGAL::Point_3<K>& q) { }

  /*!
    \ingroup PkgWeightsRefInverseDistanceWeights

    \brief computes the inverse distance weight in 2D using the points `p` and `q`,
    given a traits class `traits` with geometric objects, predicates, and constructions.
  */
  template<typename GeomTraits>
  typename GeomTraits::FT inverse_distance_weight(
    const typename GeomTraits::Point_2& p,
    const typename GeomTraits::Point_2& q,
    const GeomTraits& traits) { }

  /*!
    \ingroup PkgWeightsRefInverseDistanceWeights

    \brief computes the inverse distance weight in 3D using the points `p` and `q`,
    given a traits class `traits` with geometric objects, predicates, and constructions.
  */
  template<typename GeomTraits>
  typename GeomTraits::FT inverse_distance_weight(
    const typename GeomTraits::Point_3& p,
    const typename GeomTraits::Point_3& q,
    const GeomTraits& traits) { }

  /*!
    \ingroup PkgWeightsRefInverseDistanceWeights

    \brief computes the inverse distance weight in 2D using the points `p` and `q`,
    which are parameterized by a `Kernel` K.
  */
  template<typename K>
  typename K::FT inverse_distance_weight(
    const CGAL::Point_2<K>& p,
    const CGAL::Point_2<K>& q) { }

  /*!
    \ingroup PkgWeightsRefInverseDistanceWeights

    \brief computes the inverse distance weight in 3D using the points `p` and `q`,
    which are parameterized by a `Kernel` K.
  */
  template<typename K>
  typename K::FT inverse_distance_weight(
    const CGAL::Point_3<K>& p,
    const CGAL::Point_3<K>& q) { }

  #endif // DOXYGEN_RUNNING

  /// \cond SKIP_IN_MANUAL
  template<typename GeomTraits>
  typename GeomTraits::FT inverse_distance_weight(
    const typename GeomTraits::Point_2&,
    const typename GeomTraits::Point_2& r,
    const typename GeomTraits::Point_2&,
    const typename GeomTraits::Point_2& q,
    const GeomTraits& traits) {

    using FT = typename GeomTraits::FT;
    const FT d = internal::distance_2(traits, q, r);
    return inverse_distance_ns::weight(d);
  }

  template<typename GeomTraits>
  typename GeomTraits::FT inverse_distance_weight(
    const CGAL::Point_2<GeomTraits>& t,
    const CGAL::Point_2<GeomTraits>& r,
    const CGAL::Point_2<GeomTraits>& p,
    const CGAL::Point_2<GeomTraits>& q) {

    const GeomTraits traits;
    return inverse_distance_weight(t, r, p, q, traits);
  }

  template<typename GeomTraits>
  typename GeomTraits::FT inverse_distance_weight(
    const typename GeomTraits::Point_2& p,
    const typename GeomTraits::Point_2& q,
    const GeomTraits& traits) {

    typename GeomTraits::Point_2 stub;
    return inverse_distance_weight(stub, p, stub, q, traits);
  }

  template<typename GeomTraits>
  typename GeomTraits::FT inverse_distance_weight(
    const CGAL::Point_2<GeomTraits>& p,
    const CGAL::Point_2<GeomTraits>& q) {

    CGAL::Point_2<GeomTraits> stub;
    return inverse_distance_weight(stub, p, stub, q);
  }

  template<typename GeomTraits>
  typename GeomTraits::FT inverse_distance_weight(
    const typename GeomTraits::Point_3&,
    const typename GeomTraits::Point_3& r,
    const typename GeomTraits::Point_3&,
    const typename GeomTraits::Point_3& q,
    const GeomTraits& traits) {

    using FT = typename GeomTraits::FT;
    const FT d = internal::distance_3(traits, q, r);
    return inverse_distance_ns::weight(d);
  }

  template<typename GeomTraits>
  typename GeomTraits::FT inverse_distance_weight(
    const CGAL::Point_3<GeomTraits>& t,
    const CGAL::Point_3<GeomTraits>& r,
    const CGAL::Point_3<GeomTraits>& p,
    const CGAL::Point_3<GeomTraits>& q) {

    const GeomTraits traits;
    return inverse_distance_weight(t, r, p, q, traits);
  }

  template<typename GeomTraits>
  typename GeomTraits::FT inverse_distance_weight(
    const typename GeomTraits::Point_3& p,
    const typename GeomTraits::Point_3& q,
    const GeomTraits& traits) {

    typename GeomTraits::Point_3 stub;
    return inverse_distance_weight(stub, p, stub, q, traits);
  }

  template<typename GeomTraits>
  typename GeomTraits::FT inverse_distance_weight(
    const CGAL::Point_3<GeomTraits>& p,
    const CGAL::Point_3<GeomTraits>& q) {

    CGAL::Point_3<GeomTraits> stub;
    return inverse_distance_weight(stub, p, stub, q);
  }
  /// \endcond

} // namespace Weights
} // namespace CGAL

#endif // CGAL_INVERSE_DISTANCE_WEIGHTS_H
