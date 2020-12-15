// Copyright (c) 2019 GeometryFactory SARL (France).
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
// Author(s)     : Simon Giraudot

#ifndef CGAL_KSR_UTILS_H
#define CGAL_KSR_UTILS_H

// STL includes.
#include <set>
#include <cmath>
#include <array>
#include <string>
#include <sstream>
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

// Line discretization.
#define CGAL_KSR_SAME_VECTOR_TOLERANCE 0.99999
#define CGAL_KSR_SAME_POINT_TOLERANCE 1e-10

namespace CGAL {
namespace KSR {

// Size type.
#ifdef CGAL_KSR_USE_STD_SIZE_T_AS_SIZE_TYPE
using size_t = std::size_t;
using std::vector;
#else

using size_t = boost::uint32_t;
template<typename ValueType>
class vector {

public:
  using value_type     = ValueType;
  using Base           = std::vector<ValueType>;
  using const_iterator = typename Base::const_iterator;
  using iterator       = typename Base::iterator;

private:
  std::vector<ValueType> m_data;

public:
  vector(const KSR::size_t size = 0) :
  m_data(size)
  { }
  vector(const KSR::size_t size, const ValueType& value) :
  m_data(size, value)
  { }

  const_iterator begin() const { return m_data.begin(); }
  const_iterator end() const { return m_data.end(); }
  iterator begin() { return m_data.begin(); }
  iterator end() { return m_data.end(); }

  const KSR::size_t size() const { return static_cast<KSR::size_t>(m_data.size()); }
  const bool empty() const { return m_data.empty(); }
  void clear() { m_data.clear(); }

  void reserve(const KSR::size_t size) { m_data.reserve(std::size_t(size)); }
  void resize(const KSR::size_t size) { m_data.resize(std::size_t(size)); }

  const ValueType& operator[](const KSR::size_t idx) const { return m_data[std::size_t(idx)]; }
  ValueType& operator[](const KSR::size_t idx) { return m_data[std::size_t(idx)]; }

  void erase(const iterator it) { m_data.erase(it); }
  void insert(const iterator it, const ValueType& value) { m_data.insert(it, value); }

  const ValueType& front() const { return m_data.front(); }
  ValueType& front() { return m_data.front(); }
  const ValueType& back() const { return m_data.back(); }
  ValueType& back() { return m_data.back(); }

  void push_back(const ValueType& value) { m_data.push_back(value); }
  void swap(vector& other) { m_data.swap(other.m_data); }

  const bool operator<(const vector& other) const {
    return (this->m_data < other.m_data);
  }
};

#endif

using Idx_vector = vector<KSR::size_t>;
using Idx_vector_iterator = typename Idx_vector::iterator;

using Idx_set = std::set<KSR::size_t>;
using Idx_set_iterator = typename Idx_set::iterator;

// Use -1 as no element identifier.
inline const KSR::size_t no_element() { return KSR::size_t(-1); }

// Use -2 as special uninitialized identifier.
inline const KSR::size_t uninitialized() { return KSR::size_t(-2); }

template<typename Point_d>
const std::string to_string(const Point_d& p) {
  std::ostringstream oss;
  oss.precision(20);
  oss << p;
  return oss.str();
}

template<typename Point_d>
decltype(auto) distance(const Point_d& p, const Point_d& q) {
  using Traits = typename Kernel_traits<Point_d>::Kernel;
  using FT = typename Traits::FT;
  const FT sq_dist = CGAL::squared_distance(p, q);
  return static_cast<FT>(CGAL::sqrt(CGAL::to_double(sq_dist)));
}

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

template<typename Vector_d>
inline const Vector_d normalize(const Vector_d& v) {
  using Traits = typename Kernel_traits<Vector_d>::Kernel;
  using FT = typename Traits::FT;
  const FT dot_product = CGAL::abs(v * v);
  CGAL_assertion(dot_product != FT(0));
  return v / static_cast<FT>(CGAL::sqrt(CGAL::to_double(dot_product)));
}

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
  CGAL_assertion_msg(is_intersection_found, "ERROR: INTERSECTION IS NOT FOUND!");
  return out;
}

template<typename Segment_2>
const bool are_parallel(
  const Segment_2& seg1, const Segment_2& seg2) {

  using Traits = typename Kernel_traits<Segment_2>::Kernel;
  using FT = typename Traits::FT;

  const FT tol = tolerance<FT>();
  FT m1 = FT(100000), m2 = FT(100000);

  const FT d1 = (seg1.target().x() - seg1.source().x());
  const FT d2 = (seg2.target().x() - seg2.source().x());

  if (CGAL::abs(d1) > tol)
    m1 = (seg1.target().y() - seg1.source().y()) / d1;
  if (CGAL::abs(d2) > tol)
    m2 = (seg2.target().y() - seg2.source().y()) / d2;

  // return CGAL::parallel(seg1, seg2); // exact version

  if (CGAL::abs(m1 - m2) < tol) { // approximate version
    return true;
  }
  return false;
}

} // namespace KSR
} // namespace CGAL

#endif // CGAL_KSR_UTILS_H
