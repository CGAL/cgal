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
#include <CGAL/intersections.h>
#include <CGAL/Cartesian_converter.h>

namespace CGAL {
namespace KSR {

template<typename GeomTraits>
class Kinetic_traits_3 {

  using FT = typename GeomTraits::FT;

public:
  Kinetic_traits_3(const bool use_hybrid_mode) :
  m_use_hybrid_mode(use_hybrid_mode) { }

  template<typename ResultType, typename Type1, typename Type2>
  inline const ResultType intersection(const Type1& t1, const Type2& t2) const {

    ResultType out;
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

  template<typename Type1, typename Type2>
  decltype(auto) intersection_impl(const Type1& t1, const Type2& t2) const {

    if (!m_use_hybrid_mode) {
      return CGAL::intersection(t1, t2);
    } {
      CGAL_assertion_msg(false, "TODO: FINISH HYBRID MODE!");
      return CGAL::intersection(t1, t2);
    }
  }
};

} // namespace KSR
} // namespace CGAL

#endif // CGAL_KSR_CONVERSIONS_H
