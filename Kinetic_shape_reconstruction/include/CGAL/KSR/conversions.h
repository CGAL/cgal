// Copyright (c) 2021 GeometryFactory SARL (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Simon Giraudot, Dmitry Anisimov

#ifndef CGAL_KSR_CONVERSIONS_H
#define CGAL_KSR_CONVERSIONS_H

// #include <CGAL/license/Kinetic_shape_reconstruction.h>

// STL includes.
#include <vector>
#include <utility>

// CGAL includes.
#include <CGAL/assertions.h>
#include <CGAL/intersections.h>
#include <CGAL/Cartesian_converter.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

namespace CGAL {
namespace KSR {

template<typename GeomTraits>
class Kinetic_traits_3 {

  using IK = GeomTraits; // assume GT is inexact here

  using IT      = typename GeomTraits::FT;
  using Point_2 = typename GeomTraits::Point_2;
  using Line_2  = typename GeomTraits::Line_2;

  // TODO: This is very naive way of doing this. We should compare IT and ET
  // and, in case they are the same or IT is already exact, abort conversion!
  using EK = CGAL::Exact_predicates_exact_constructions_kernel;
  using ET = typename EK::FT;

  using IK_to_EK = CGAL::Cartesian_converter<IK, EK>;
  using EK_to_IK = CGAL::Cartesian_converter<EK, IK>;

public:
  Kinetic_traits_3() :
  m_use_hybrid_mode(false) { }

  inline const Point_2 intersection(const Line_2& t1, const Line_2& t2) const {

    Point_2 out;
    const bool is_intersection_found = intersection(t1, t2, out);
    CGAL_assertion(is_intersection_found);
    return out;
  }

  template<typename Type1, typename Type2, typename ResultType>
  inline bool intersection(
    const Type1& t1, const Type2& t2, ResultType& result) const {

    const auto inter = intersection_impl(t1, t2);
    if (!inter) return false;
    if (const ResultType* typed_inter = boost::get<ResultType>(&*inter)) {
      result = *typed_inter;
      return true;
    }
    return false;
  }

private:
  const bool m_use_hybrid_mode;
  IK_to_EK m_inexact_to_exact;
  EK_to_IK m_exact_to_inexact;

  template<typename Type1, typename Type2>
  decltype(auto) intersection_impl(const Type1& t1, const Type2& t2) const {

    // if (!m_use_hybrid_mode) {
      return CGAL::intersection(t1, t2);
    // } else { // convert to exact

    //   // TODO: It does not compile with EPECK as input kernel.
    //   const auto exact_t1 = m_inexact_to_exact(t1);
    //   const auto exact_t2 = m_inexact_to_exact(t2);
    //   const auto exact_result = CGAL::intersection(exact_t1, exact_t2);
    //   return m_exact_to_inexact(exact_result);
    // }
  }
};

} // namespace KSR
} // namespace CGAL

#endif // CGAL_KSR_CONVERSIONS_H
