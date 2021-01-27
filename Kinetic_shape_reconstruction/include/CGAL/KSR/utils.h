// Copyright (c) 2020 GeometryFactory SARL (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
//
// Author(s)     : Dmitry Anisimov, Simon Giraudot

#ifndef CGAL_KSR_UTILS_H
#define CGAL_KSR_UTILS_H

// #include <CGAL/license/Kinetic_shape_reconstruction.h>

// STL includes.
#include <set>
#include <cmath>
#include <array>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <deque>
#include <queue>
#include <map>

// CGAL includes.
#include <CGAL/Bbox_3.h>
#include <CGAL/Iterator_range.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/centroid.h>

// Boost includes.
#include <boost/function_output_iterator.hpp>

namespace CGAL {
namespace KSR {

// Use -1 as no element identifier.
inline const std::size_t no_element() { return std::size_t(-1); }

// Use -2 as special uninitialized identifier.
inline const std::size_t uninitialized() { return std::size_t(-2); }

// Convert point to string.
template<typename Point_d>
const std::string to_string(const Point_d& p) {
  std::ostringstream oss;
  oss.precision(20);
  oss << p;
  return oss.str();
}

// Distance between two points.
template<typename Point_d>
decltype(auto) distance(const Point_d& p, const Point_d& q) {
  using Traits = typename Kernel_traits<Point_d>::Kernel;
  using FT = typename Traits::FT;
  const FT sq_dist = CGAL::squared_distance(p, q);
  return static_cast<FT>(CGAL::sqrt(CGAL::to_double(sq_dist)));
}

// Project 3D point onto 2D plane.
template<typename Point_3>
typename Kernel_traits<Point_3>::Kernel::Point_2
point_2_from_point_3(const Point_3& point_3) {
  return typename Kernel_traits<Point_3>::Kernel::Point_2(
    point_3.x(), point_3.y());
}

// Tolerance.
template<typename FT>
static FT tolerance() {
  return FT(1) / FT(100000);
}

template<typename FT>
static FT point_tolerance() {
  return tolerance<FT>();
}

template<typename FT>
static FT vector_tolerance() {
  return FT(99999) / FT(100000);
}

// Normalize vector.
template<typename Vector_d>
inline const Vector_d normalize(const Vector_d& v) {
  using Traits = typename Kernel_traits<Vector_d>::Kernel;
  using FT = typename Traits::FT;
  const FT dot_product = CGAL::abs(v * v);
  CGAL_assertion(dot_product != FT(0));
  return v / static_cast<FT>(CGAL::sqrt(CGAL::to_double(dot_product)));
}

// Intersections.
template<typename Type1, typename Type2, typename ResultType>
inline const bool intersection(
  const Type1& t1, const Type2& t2, ResultType& result) {

  const auto inter = intersection(t1, t2);
  if (!inter) return false;
  if (const ResultType* typed_inter = boost::get<ResultType>(&*inter)) {
    result = *typed_inter;
    return true;
  }
  return false;
}

template<typename ResultType, typename Type1, typename Type2>
inline const ResultType intersection(const Type1& t1, const Type2& t2) {

  ResultType out;
  const bool is_intersection_found = intersection(t1, t2, out);
  CGAL_assertion(is_intersection_found);
  return out;
}

// Predicates.
template<typename Segment_2>
const bool are_parallel(
  const Segment_2& seg1, const Segment_2& seg2) {

  using Traits = typename Kernel_traits<Segment_2>::Kernel;
  using FT = typename Traits::FT;

  const FT tol = tolerance<FT>();
  FT m1 = FT(100000), m2 = FT(100000);

  const FT d1 = (seg1.target().x() - seg1.source().x());
  const FT d2 = (seg2.target().x() - seg2.source().x());

  if (CGAL::abs(d1) > tol) {
    CGAL_assertion(d1 != FT(0));
    m1 = (seg1.target().y() - seg1.source().y()) / d1;
  }
  if (CGAL::abs(d2) > tol) {
    CGAL_assertion(d2 != FT(0));
    m2 = (seg2.target().y() - seg2.source().y()) / d2;
  }

  // return CGAL::parallel(seg1, seg2); // exact version
  if (CGAL::abs(m1 - m2) < tol) { // approximate version
    return true;
  }
  return false;
}

// Helpers.
template<typename IVertex>
class Indexer {
public:
  const std::size_t operator()(const IVertex& ivertex) {
    const auto pair = m_indices.insert(
      std::make_pair(ivertex, m_indices.size()));
    const auto& item = pair.first;
    const std::size_t idx = item->second;
    return idx;
  }
  void clear() { m_indices.clear(); }

private:
  std::map<IVertex, std::size_t> m_indices;
};

} // namespace KSR
} // namespace CGAL

#endif // CGAL_KSR_UTILS_H
