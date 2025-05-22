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
// Author(s)     : Simon Giraudot, Dmitry Anisimov

#ifndef CGAL_KSP_UTILS_H
#define CGAL_KSP_UTILS_H

#include <CGAL/license/Kinetic_space_partition.h>

// STL includes.
#include <set>
#include <cmath>
#include <array>
#include <string>
#include <sstream>
#include <functional>
#include <fstream>
#include <vector>
#include <deque>
#include <queue>
#include <map>

// CGAL includes.
#include <CGAL/Iterator_range.h>
#include <CGAL/number_utils.h>
#include <CGAL/assertions.h>

#include <CGAL/Cartesian_converter.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

// Boost includes.
#include <boost/iterator/function_output_iterator.hpp>

namespace CGAL {
namespace KSP {
namespace internal {

#ifdef DOXYGEN_RUNNING
#else

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
  return static_cast<FT>(CGAL::approximate_sqrt(sq_dist));
}

// Normalize vector.
template<typename Vector_d>
inline const Vector_d normalize(const Vector_d& v) {
  using Traits = typename Kernel_traits<Vector_d>::Kernel;
  using FT = typename Traits::FT;
  const FT dot_product = CGAL::abs(v * v);
  //CGAL_assertion(dot_product != FT(0));
  return v / static_cast<FT>(CGAL::approximate_sqrt(dot_product));
}

// Intersections. Used only in the 2D version.
// For the 3D version, see conversions.h!
template<typename Type1, typename Type2, typename ResultType>
inline bool intersection(
  const Type1& t1, const Type2& t2, ResultType& result) {

  const auto inter = intersection(t1, t2);
  if (!inter) return false;

  if (CGAL::assign(result, inter))
    return true;

  return false;
}

template<typename ResultType, typename Type1, typename Type2>
inline const ResultType intersection(const Type1& t1, const Type2& t2) {

  ResultType out;
  CGAL_assertion_code(const bool is_intersection_found =) intersection(t1, t2, out);
  CGAL_assertion(is_intersection_found);
  return out;
}

// Angles.

// Converts radians to degrees.
template<typename FT>
const FT degrees_2(const FT angle_rad) {
  return angle_rad * FT(180) / static_cast<FT>(CGAL_PI);
}

// Computes an angle in degrees between two directions.
template<typename Direction_2>
const typename Kernel_traits<Direction_2>::Kernel::FT
compute_angle_2(const Direction_2& dir1, const Direction_2& dir2) {

  using Traits = typename Kernel_traits<Direction_2>::Kernel;
  using FT = typename Traits::FT;

  const auto v1 = dir2.to_vector();
  const auto v2 = -dir1.to_vector();

  const FT det = CGAL::determinant(v1, v2);
  const FT dot = CGAL::scalar_product(v1, v2);
  const FT angle_rad = static_cast<FT>(
    std::atan2(CGAL::to_double(det), CGAL::to_double(dot)));
  const FT angle_deg = degrees_2(angle_rad);
  return angle_deg;
}

// Converts an angle in degrees from the range [-180, 180]
// into the mod 90 angle.
template<typename FT>
const FT convert_angle_2(const FT angle_2) {

  FT angle = angle_2;
  if (angle > FT(90)) angle = FT(180) - angle;
  else if (angle < -FT(90)) angle = FT(180) + angle;
  return angle;
}

// Computes a positive angle in degrees that
// is always in the range [0, 90].
template<typename Direction_2>
const typename Kernel_traits<Direction_2>::Kernel::FT
angle_2(const Direction_2& dir1, const Direction_2& dir2) {
  const auto angle_2 = compute_angle_2(dir1, dir2);
  return CGAL::abs(convert_angle_2(angle_2));
}

#endif

} // namespace internal
} // namespace KSP
} // namespace CGAL

#endif // CGAL_KSP_UTILS_H
