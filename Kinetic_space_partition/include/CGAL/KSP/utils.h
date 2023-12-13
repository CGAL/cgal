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
#include <CGAL/Bbox_3.h>
#include <CGAL/centroid.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Iterator_range.h>
#include <CGAL/convex_hull_2.h>
#include <CGAL/number_utils.h>
#include <CGAL/assertions.h>

#include <CGAL/Cartesian_converter.h>
#include <CGAL/linear_least_squares_fitting_2.h>
#include <CGAL/linear_least_squares_fitting_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_2_algorithms.h>

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
  return static_cast<FT>(CGAL::sqrt(CGAL::to_double(sq_dist)));
}

// Project 3D point onto 2D plane.
template<typename Point_3>
typename Kernel_traits<Point_3>::Kernel::Point_2
point_2_from_point_3(const Point_3& point_3) {
  return typename Kernel_traits<Point_3>::Kernel::Point_2(
    point_3.x(), point_3.y());
}

// Get 3D point from a 2D point.
template<typename Point_2>
typename Kernel_traits<Point_2>::Kernel::Point_3
point_3_from_point_2(const Point_2& point_2) {
  return typename Kernel_traits<Point_2>::Kernel::Point_3(
    point_2.x(), point_2.y(), typename Kernel_traits<Point_2>::Kernel::FT(0));
}

// Normalize vector.
template<typename Vector_d>
inline const Vector_d normalize(const Vector_d& v) {
  using Traits = typename Kernel_traits<Vector_d>::Kernel;
  using FT = typename Traits::FT;
  const FT dot_product = CGAL::abs(v * v);
  //CGAL_assertion(dot_product != FT(0));
  return v / static_cast<FT>(CGAL::sqrt(CGAL::to_double(dot_product)));
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
  const bool is_intersection_found = intersection(t1, t2, out);
  CGAL_assertion(is_intersection_found);
  return out;
}

// Get boundary points from a set of points.
template<typename Point_2, typename Line_2>
void boundary_points_on_line_2(
  const std::vector<Point_2>& input_range,
  const std::vector<std::size_t>& indices,
  const Line_2& line, Point_2& p, Point_2& q) {

  using Traits = typename Kernel_traits<Point_2>::Kernel;
  using FT = typename Traits::FT;
  using Vector_2 = typename Traits::Vector_2;

  FT min_proj_value = (std::numeric_limits<FT>::max)();
  FT max_proj_value = -min_proj_value;

  const auto ref_vector = line.to_vector();
  const auto& ref_point = input_range[indices.front()];

  for (const std::size_t index : indices) {
    const auto& query = input_range[index];
    const auto point = line.projection(query);
    const Vector_2 curr_vector(ref_point, point);
    const FT value = CGAL::scalar_product(curr_vector, ref_vector);

    if (value < min_proj_value) {
      min_proj_value = value;
      p = point;
    }
    if (value > max_proj_value) {
      max_proj_value = value;
      q = point;
    }
  }
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

// Classes.
template<typename IVertex>
class Indexer {
public:
  std::size_t operator()(const IVertex& ivertex) {
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

template<
  typename GeomTraits,
  typename InputRange,
  typename NeighborQuery>
class Estimate_normals_2 {

public:
  using Traits = GeomTraits;
  using Input_range = InputRange;
  using Neighbor_query = NeighborQuery;

  using Kernel = Traits;
  using FT = typename Kernel::FT;
  using Vector_2 = typename Kernel::Vector_2;
  using Line_2 = typename Kernel::Line_2;

  using Indices = std::vector<std::size_t>;

  using IK = CGAL::Exact_predicates_inexact_constructions_kernel;
  using IPoint_2 = typename IK::Point_2;
  using ILine_2 = typename IK::Line_2;
  using Converter = CGAL::Cartesian_converter<Kernel, IK>;

  Estimate_normals_2(
    const Input_range& input_range,
    const Neighbor_query& neighbor_query) :
    m_input_range(input_range),
    m_neighbor_query(neighbor_query) {

    CGAL_precondition(input_range.size() > 0);
  }

  void get_normals(std::vector<Vector_2>& normals) const {

    normals.clear();
    normals.reserve(m_input_range.size());

    Indices neighbors;
    for (std::size_t i = 0; i < m_input_range.size(); ++i) {
      neighbors.clear();
      m_neighbor_query(i, neighbors);
      const auto line = fit_line(neighbors);
      auto normal = line.to_vector();
      normal = normal.perpendicular(CGAL::COUNTERCLOCKWISE);
      normal = normalize(normal);
      normals.push_back(normal);
    }
    CGAL_assertion(normals.size() == m_input_range.size());
  }

private:
  const Input_range& m_input_range;
  const Neighbor_query& m_neighbor_query;
  const Converter m_converter;

  const Line_2 fit_line(const Indices& indices) const {
    CGAL_assertion(indices.size() > 0);

    std::vector<IPoint_2> points;
    points.reserve(indices.size());
    for (const std::size_t index : indices) {
      const auto& point = get(m_neighbor_query.point_map(), index);
      points.push_back(m_converter(point));
    }
    CGAL_assertion(points.size() == indices.size());

    ILine_2 fitted_line;
    IPoint_2 fitted_centroid;
    CGAL::linear_least_squares_fitting_2(
      points.begin(), points.end(),
      fitted_line, fitted_centroid,
      CGAL::Dimension_tag<0>());

    const Line_2 line(
      static_cast<FT>(fitted_line.a()),
      static_cast<FT>(fitted_line.b()),
      static_cast<FT>(fitted_line.c()));
    return line;
  }
};

#endif

} // namespace internal
} // namespace KSP
} // namespace CGAL

#endif // CGAL_KSP_UTILS_H
