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

#ifndef CGAL_BARYCENTRIC_INTERNAL_COTANGENT_WEIGHTS_2_H
#define CGAL_BARYCENTRIC_INTERNAL_COTANGENT_WEIGHTS_2_H

#include <CGAL/license/Barycentric_coordinates_2.h>

namespace CGAL {
namespace Barycentric_coordinates {
namespace internal {

  namespace cotangent_ns {

    template<typename FT>
    const FT weight(
      const FT cot_beta, const FT cot_gamma) {

      return FT(2) * (cot_beta + cot_gamma);
    }
  }

  // Computes cotanget between two 2D vectors.
  template<typename GeomTraits>
  const typename GeomTraits::FT cotangent_2(
    const GeomTraits& traits,
    const typename GeomTraits::Point_2& p,
    const typename GeomTraits::Point_2& q,
    const typename GeomTraits::Point_2& r) {

    using FT = typename GeomTraits::FT;
    const auto dot_product_2 =
      traits.compute_scalar_product_2_object();
    const auto cross_product_2 =
      traits.compute_determinant_2_object();
    const auto construct_vector_2 =
      traits.construct_vector_2_object();

    const auto v1 = construct_vector_2(q, r);
    const auto v2 = construct_vector_2(q, p);

    const FT dot = dot_product_2(v1, v2);
    const FT cross = cross_product_2(v1, v2);

    const FT length = CGAL::abs(cross);
    CGAL_assertion(length != FT(0));
    if (length != FT(0)) return dot / length;
    else return FT(0); // undefined
  }

  template<typename GeomTraits>
  const typename GeomTraits::FT cotangent_weight(
    const typename GeomTraits::Point_2& t,
    const typename GeomTraits::Point_2& r,
    const typename GeomTraits::Point_2& p,
    const typename GeomTraits::Point_2& q,
    const GeomTraits& traits) {

    using FT = typename GeomTraits::FT;
    const FT cot_beta  = cotangent_2(traits, q, t, r);
    const FT cot_gamma = cotangent_2(traits, r, p, q);
    return cotangent_ns::weight(cot_beta, cot_gamma);
  }

} // namespace internal
} // namespace Barycentric_coordinates
} // namespace CGAL

#endif // CGAL_BARYCENTRIC_INTERNAL_COTANGENT_WEIGHTS_2_H
