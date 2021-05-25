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

#ifndef CGAL_THREE_POINT_FAMILY_WEIGHTS_H
#define CGAL_THREE_POINT_FAMILY_WEIGHTS_H

// #include <CGAL/license/Weights.h>

// Internal includes.
#include <CGAL/Weights/internal/utils.h>

namespace CGAL {
namespace Weights {

  /// \cond SKIP_IN_MANUAL
  namespace three_point_family_ns {

    template<typename GeomTraits>
    typename GeomTraits::FT weight(
      const GeomTraits& traits,
      const typename GeomTraits::FT d1,
      const typename GeomTraits::FT d2,
      const typename GeomTraits::FT d3,
      const typename GeomTraits::FT A1,
      const typename GeomTraits::FT A2,
      const typename GeomTraits::FT B,
      const typename GeomTraits::FT p) {

      using FT = typename GeomTraits::FT;
      FT w = FT(0);
      CGAL_precondition(A1 != FT(0) && A2 != FT(0));
      const FT prod = A1 * A2;
      if (prod != FT(0)) {
        const FT inv = FT(1) / prod;
        FT r1 = d1;
        FT r2 = d2;
        FT r3 = d3;
        if (p != FT(1)) {
          r1 = internal::power(traits, d1, p);
          r2 = internal::power(traits, d2, p);
          r3 = internal::power(traits, d3, p);
        }
        w = (r3 * A1 - r2 * B + r1 * A2) * inv;
      }
      return w;
    }
  }
  /// \endcond

  #if defined(DOXYGEN_RUNNING)

  /*!
    \ingroup PkgWeightsRefThreePointFamilyWeights

    \brief computes the three-point family weight in 2D at `q` using the points `p0`, `p1`,
    and `p2` and the power parameter `a`, given a traits class `traits` with geometric objects,
    predicates, and constructions.
  */
  template<typename GeomTraits>
  typename GeomTraits::FT three_point_family_weight(
    const typename GeomTraits::Point_2& p0,
    const typename GeomTraits::Point_2& p1,
    const typename GeomTraits::Point_2& p2,
    const typename GeomTraits::Point_2& q,
    const typename GeomTraits::FT a,
    const GeomTraits& traits) { }

  /*!
    \ingroup PkgWeightsRefThreePointFamilyWeights

    \brief computes the three-point family weight in 2D at `q` using the points `p0`, `p1`,
    and `p2`, which are parameterized by a `Kernel` K, and the power parameter `a`, which
    can be omitted.
  */
  template<typename K>
  typename K::FT three_point_family_weight(
    const CGAL::Point_2<K>& p0,
    const CGAL::Point_2<K>& p1,
    const CGAL::Point_2<K>& p2,
    const CGAL::Point_2<K>& q,
    const typename K::FT a = typename K::FT(1)) { }

  #endif // DOXYGEN_RUNNING

  /// \cond SKIP_IN_MANUAL
  template<typename GeomTraits>
  typename GeomTraits::FT three_point_family_weight(
    const typename GeomTraits::Point_2& t,
    const typename GeomTraits::Point_2& r,
    const typename GeomTraits::Point_2& p,
    const typename GeomTraits::Point_2& q,
    const typename GeomTraits::FT a,
    const GeomTraits& traits) {

    using FT = typename GeomTraits::FT;
    const FT d1 = internal::distance_2(traits, q, t);
    const FT d2 = internal::distance_2(traits, q, r);
    const FT d3 = internal::distance_2(traits, q, p);

    const FT A1 = internal::area_2(traits, r, q, t);
    const FT A2 = internal::area_2(traits, p, q, r);
    const FT B  = internal::area_2(traits, p, q, t);

    return three_point_family_ns::weight(
      traits, d1, d2, d3, A1, A2, B, a);
  }

  template<typename GeomTraits>
  typename GeomTraits::FT three_point_family_weight(
    const CGAL::Point_2<GeomTraits>& t,
    const CGAL::Point_2<GeomTraits>& r,
    const CGAL::Point_2<GeomTraits>& p,
    const CGAL::Point_2<GeomTraits>& q,
    const typename GeomTraits::FT a =
    typename GeomTraits::FT(1)) {

    const GeomTraits traits;
    return three_point_family_weight(t, r, p, q, a, traits);
  }

  namespace internal {

  // Example of flattening:

  // 3D configuration.
  // const Point_3 p0(0, 1, 1);
  // const Point_3 p1(2, 0, 1);
  // const Point_3 p2(7, 1, 1);
  // const Point_3 q0(3, 1, 1);

  // Choose a type of the weight:
  // e.g. 0 - Wachspress (WP) weight.
  // const FT wp = FT(0);

  // Compute WP weights for q1, which is not on the plane [p0, p1, p2].

  // Point_3 q1(3, 1, 2);
  // std::cout << "3D wachspress (WP, q1): ";
  // std::cout << CGAL::Weights::three_point_family_weight(p0, p1, p2, q1, wp) << std::endl;

  // Converge q1 towards q0 that is we flatten the configuration.
  // We also compare the result with the authalic weight.

  // std::cout << "Converge q1 to q0: " << std::endl;
  // for (FT x = FT(0); x <= FT(1); x += step) {
  //   std::cout << "3D wachspress/authalic: ";
  //   q1 = Point_3(3, 1, FT(2) - x);
  //   std::cout << CGAL::Weights::three_point_family_weight(p0, p1, p2, q1, wp) << "/";
  //   std::cout << CGAL::Weights::authalic_weight(p0, p1, p2, q1) << std::endl;
  // }

  template<typename GeomTraits>
  typename GeomTraits::FT three_point_family_weight(
    const typename GeomTraits::Point_3& t,
    const typename GeomTraits::Point_3& r,
    const typename GeomTraits::Point_3& p,
    const typename GeomTraits::Point_3& q,
    const typename GeomTraits::FT a,
    const GeomTraits& traits) {

    using Point_2 = typename GeomTraits::Point_2;
    Point_2 tf, rf, pf, qf;
    internal::flatten(
      traits,
      t,  r,  p,  q,
      tf, rf, pf, qf);
    return CGAL::Weights::
      three_point_family_weight(tf, rf, pf, qf, a, traits);
  }

  template<typename GeomTraits>
  typename GeomTraits::FT three_point_family_weight(
    const CGAL::Point_3<GeomTraits>& t,
    const CGAL::Point_3<GeomTraits>& r,
    const CGAL::Point_3<GeomTraits>& p,
    const CGAL::Point_3<GeomTraits>& q,
    const typename GeomTraits::FT a =
    typename GeomTraits::FT(1)) {

    const GeomTraits traits;
    return three_point_family_weight(t, r, p, q, a, traits);
  }

  } // namespace internal

  /// \endcond

} // namespace Weights
} // namespace CGAL

#endif // CGAL_THREE_POINT_FAMILY_WEIGHTS_H
