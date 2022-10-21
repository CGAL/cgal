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

#ifndef CGAL_UNIFORM_REGION_WEIGHTS_H
#define CGAL_UNIFORM_REGION_WEIGHTS_H

// Internal includes.
#include <CGAL/Weights/internal/utils.h>

namespace CGAL {
namespace Weights {

  #if defined(DOXYGEN_RUNNING)

  /*!
    \ingroup PkgWeightsRefUniformRegionWeights

    \brief this function always returns 1, given three points in 2D and a traits class
    with geometric objects, predicates, and constructions.
  */
  template<typename GeomTraits>
  typename GeomTraits::FT uniform_area(
    const typename GeomTraits::Point_2&,
    const typename GeomTraits::Point_2&,
    const typename GeomTraits::Point_2&,
    const GeomTraits&) { }

  /*!
    \ingroup PkgWeightsRefUniformRegionWeights

    \brief this function always returns 1, given three points in 3D and a traits class
    with geometric objects, predicates, and constructions.
  */
  template<typename GeomTraits>
  typename GeomTraits::FT uniform_area(
    const typename GeomTraits::Point_3&,
    const typename GeomTraits::Point_3&,
    const typename GeomTraits::Point_3&,
    const GeomTraits&) { }

  /*!
    \ingroup PkgWeightsRefUniformRegionWeights

    \brief this function always returns 1, given three points in 2D which are
    parameterized by a `Kernel` K.
  */
  template<typename K>
  typename K::FT uniform_area(
    const CGAL::Point_2<K>&,
    const CGAL::Point_2<K>&,
    const CGAL::Point_2<K>&) { }

  /*!
    \ingroup PkgWeightsRefUniformRegionWeights

    \brief this function always returns 1, given three points in 3D which are
    parameterized by a `Kernel` K.
  */
  template<typename K>
  typename K::FT uniform_area(
    const CGAL::Point_3<K>&,
    const CGAL::Point_3<K>&,
    const CGAL::Point_3<K>&) { }

  #endif // DOXYGEN_RUNNING

  /// \cond SKIP_IN_MANUAL
  template<typename GeomTraits>
  typename GeomTraits::FT uniform_area(
    const typename GeomTraits::Point_2&,
    const typename GeomTraits::Point_2&,
    const typename GeomTraits::Point_2&,
    const GeomTraits&) {

    using FT = typename GeomTraits::FT;
    return FT(1);
  }

  template<typename GeomTraits>
  typename GeomTraits::FT uniform_area(
    const CGAL::Point_2<GeomTraits>& p,
    const CGAL::Point_2<GeomTraits>& q,
    const CGAL::Point_2<GeomTraits>& r) {

    const GeomTraits traits;
    return uniform_area(p, q, r, traits);
  }

  template<typename GeomTraits>
  typename GeomTraits::FT uniform_area(
    const typename GeomTraits::Point_3&,
    const typename GeomTraits::Point_3&,
    const typename GeomTraits::Point_3&,
    const GeomTraits&) {

    using FT = typename GeomTraits::FT;
    return FT(1);
  }

  template<typename GeomTraits>
  typename GeomTraits::FT uniform_area(
    const CGAL::Point_3<GeomTraits>& p,
    const CGAL::Point_3<GeomTraits>& q,
    const CGAL::Point_3<GeomTraits>& r) {

    const GeomTraits traits;
    return uniform_area(p, q, r, traits);
  }
  /// \endcond

} // namespace Weights
} // namespace CGAL

#endif // CGAL_UNIFORM_REGION_WEIGHTS_H
