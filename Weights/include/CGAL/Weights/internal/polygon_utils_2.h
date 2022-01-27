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

#ifndef CGAL_WEIGHTS_INTERNAL_POLYGON_UTILS_2_H
#define CGAL_WEIGHTS_INTERNAL_POLYGON_UTILS_2_H

// STL includes.
#include <cmath>
#include <string>
#include <memory>
#include <vector>
#include <utility>
#include <iterator>
#include <iostream>

// CGAL includes.
#include <CGAL/utils.h>
#include <CGAL/assertions.h>
#include <CGAL/number_utils.h>
#include <CGAL/property_map.h>
#include <CGAL/Polygon_2_algorithms.h>

namespace CGAL {
namespace Weights {
namespace internal {

  enum class Edge_case {

    EXTERIOR = 0, // exterior part of the polygon
    BOUNDARY = 1, // boundary part of the polygon
    INTERIOR = 2  // interior part of the polygon
  };

  // VertexRange type enum.
  enum class Polygon_type {

    // Concave polygon = non-convex polygon.
    CONCAVE = 0,

    // This is a convex polygon with collinear vertices.
    WEAKLY_CONVEX = 1,

    // This is a convex polygon without collinear vertices.
    STRICTLY_CONVEX = 2
  };

  // This function is taken from the Polygon_2_algorithms.h header.
  // But it is updated to support property maps.
  template<
  class Point_2,
  class Orientation_2,
  class CompareX_2>
  int which_side_in_slab_2(
    const Point_2& query, const Point_2& low, const Point_2& high,
    const Orientation_2& orientation_2, const CompareX_2& compare_x_2) {

    const auto low_x_comp_res = compare_x_2(query, low);
    const auto high_x_comp_res = compare_x_2(query, high);
    if (low_x_comp_res == CGAL::SMALLER) {
      if (high_x_comp_res == CGAL::SMALLER) {
        return -1;
      }
    } else {
      switch (high_x_comp_res) {
        case CGAL::LARGER: return 1;
        case CGAL::SMALLER: break;
        case CGAL::EQUAL: return (low_x_comp_res == CGAL::EQUAL) ? 0 : 1;
      }
    }
    switch (orientation_2(low, query, high)) {
      case CGAL::LEFT_TURN:  return  1;
      case CGAL::RIGHT_TURN: return -1;
      default: return 0;
    }
  }

  // This function is taken from the Polygon_2_algorithms.h header.
  // But it is updated to support property maps.
  template<
  typename VertexRange,
  typename GeomTraits,
  typename PointMap>
  Edge_case bounded_side_2(
    const VertexRange& polygon, const typename GeomTraits::Point_2& query,
    const GeomTraits& traits, const PointMap point_map) {

    const auto first = polygon.begin();
    const auto last  = polygon.end();

    auto curr = first;
    if (curr == last) {
      return Edge_case::EXTERIOR;
    }

    auto next = curr; ++next;
    if (next == last) {
      return Edge_case::EXTERIOR;
    }

    const auto compare_x_2 = traits.compare_x_2_object();
    const auto compare_y_2 = traits.compare_y_2_object();
    const auto orientation_2 = traits.orientation_2_object();

    bool is_inside = false;
    auto curr_y_comp_res = compare_y_2(get(point_map, *curr), query);

    // Check if the segment (curr, next) intersects
    // the ray { (t, query.y()) | t >= query.x() }.
    do {
      const auto& currp = get(point_map, *curr);
      const auto& nextp = get(point_map, *next);

      auto next_y_comp_res = compare_y_2(nextp, query);
      switch (curr_y_comp_res) {
        case CGAL::SMALLER:
          switch (next_y_comp_res) {
            case CGAL::SMALLER:
              break;
            case CGAL::EQUAL:
              switch (compare_x_2(query, nextp)) {
                case CGAL::SMALLER: is_inside = !is_inside; break;
                case CGAL::EQUAL:   return Edge_case::BOUNDARY;
                case CGAL::LARGER:  break;
              }
              break;
            case CGAL::LARGER:
              switch (which_side_in_slab_2(
                query, currp, nextp, orientation_2, compare_x_2)) {
                case -1: is_inside = !is_inside; break;
                case  0: return Edge_case::BOUNDARY;
              }
              break;
          }
          break;
        case CGAL::EQUAL:
          switch (next_y_comp_res) {
            case CGAL::SMALLER:
              switch (compare_x_2(query, currp)) {
                case CGAL::SMALLER: is_inside = !is_inside; break;
                case CGAL::EQUAL:   return Edge_case::BOUNDARY;
                case CGAL::LARGER:  break;
              }
              break;
            case CGAL::EQUAL:
              switch (compare_x_2(query, currp)) {
                case CGAL::SMALLER:
                  if (compare_x_2(query, nextp) != CGAL::SMALLER) {
                    return Edge_case::BOUNDARY;
                  }
                  break;
                case CGAL::EQUAL: return Edge_case::BOUNDARY;
                case CGAL::LARGER:
                  if (compare_x_2(query, nextp) != CGAL::LARGER) {
                    return Edge_case::BOUNDARY;
                  }
                  break;
              }
              break;
            case CGAL::LARGER:
              if (compare_x_2(query, currp) == CGAL::EQUAL) {
                return Edge_case::BOUNDARY;
              }
              break;
          }
          break;
        case CGAL::LARGER:
          switch (next_y_comp_res) {
            case CGAL::SMALLER:
              switch (which_side_in_slab_2(
                query, nextp, currp, orientation_2, compare_x_2)) {
                case -1: is_inside = !is_inside; break;
                case  0: return Edge_case::BOUNDARY;
              }
              break;
            case CGAL::EQUAL:
              if (compare_x_2(query, nextp) == CGAL::EQUAL) {
                return Edge_case::BOUNDARY;
              }
              break;
            case CGAL::LARGER:
              break;
          }
          break;
      }

      curr = next;
      curr_y_comp_res = next_y_comp_res;
      ++next;
      if (next == last) {
        next = first;
      }
    } while (curr != first);
    return is_inside ? Edge_case::INTERIOR : Edge_case::EXTERIOR;
  }

  // This function is taken from the Polygon_2_algorithms.h header.
  // But it is updated to support property maps.
  template<
  typename VertexRange,
  typename GeomTraits,
  typename PointMap>
  bool is_convex_2(
    const VertexRange& polygon, const GeomTraits traits, const PointMap point_map) {

    auto first = polygon.begin();
    const auto last  = polygon.end();

    auto prev = first;
    if (prev == last) {
      return true;
    }

    auto curr = prev; ++curr;
    if (curr == last) {
      return true;
    }

    auto next = curr; ++next;
    if (next == last) {
      return true;
    }

    const auto equal_2 = traits.equal_2_object();
    while (equal_2(get(point_map, *prev), get(point_map, *curr))) {
      curr = next; ++next;
      if (next == last) {
        return true;
      }
    }

    const auto less_xy_2 = traits.less_xy_2_object();
    const auto orientation_2 = traits.orientation_2_object();

    bool has_clockwise_triplets = false;
    bool has_counterclockwise_triplets = false;
    bool order = less_xy_2(
      get(point_map, *prev), get(point_map, *curr));
    int num_order_changes = 0;

    do {
    switch_orient:
      switch (orientation_2(
        get(point_map, *prev), get(point_map, *curr), get(point_map, *next))) {

        case CGAL::CLOCKWISE:
          has_clockwise_triplets = true;
          break;
        case CGAL::COUNTERCLOCKWISE:
          has_counterclockwise_triplets = true;
          break;
        case CGAL::ZERO: {
          if (equal_2(
            get(point_map, *curr),
            get(point_map, *next))) {

            if (next == first) {
              first = curr;
            }
            ++next;
            if (next == last) {
              next = first;
            }
            goto switch_orient;
          }
          break;
        }
      }

      const bool new_order = less_xy_2(
        get(point_map, *curr), get(point_map, *next));

      if (order != new_order) {
        num_order_changes++;
      }

      if (num_order_changes > 2) {
        return false;
      }

      if (has_clockwise_triplets && has_counterclockwise_triplets) {
        return false;
      }

      prev = curr;
      curr = next;
      ++next;
      if (next == last) {
        next = first;
      }
      order = new_order;
    } while (prev != first);
    return true;
  }

  // This function is taken from the Polygon_2_algorithms.h header.
  // But it is updated to support property maps.
  template<
  typename VertexRange,
  typename GeomTraits,
  typename PointMap>
  bool is_simple_2(
    const VertexRange& polygon, const GeomTraits traits, const PointMap point_map) {

    const auto first = polygon.begin();
    const auto last = polygon.end();
    if (first == last) {
      return true;
    }

    std::vector<typename GeomTraits::Point_2> poly;
    poly.reserve(polygon.size());
    for (const auto& vertex : polygon) {
      poly.push_back(get(point_map, vertex));
    }
    return CGAL::is_simple_2(poly.begin(), poly.end(), traits);
  }

  template<
  typename VertexRange,
  typename GeomTraits,
  typename PointMap>
  Polygon_type polygon_type_2(
    const VertexRange& polygon, const GeomTraits traits, const PointMap point_map) {

    const auto collinear_2 =
      traits.collinear_2_object();
    CGAL_precondition(polygon.size() >= 3);

    // First, test the polygon on convexity.
    if (is_convex_2(polygon, traits, point_map)) {

      // Test all the consequent triplets of polygon vertices on collinearity.
      // In case we find at least one, return WEAKLY_CONVEX polygon.
      const std::size_t n = polygon.size();
      for (std::size_t i = 0; i < n; ++i) {
        const auto& p1 = get(point_map, *(polygon.begin() + i));

        const std::size_t im = (i + n - 1) % n;
        const std::size_t ip = (i + 1) % n;

        const auto& p0 = get(point_map, *(polygon.begin() + im));
        const auto& p2 = get(point_map, *(polygon.begin() + ip));

        if (collinear_2(p0, p1, p2)) {
          return Polygon_type::WEAKLY_CONVEX;
        }
      }
      // Otherwise, return STRICTLY_CONVEX polygon.
      return Polygon_type::STRICTLY_CONVEX;
    }
    // Otherwise, return CONCAVE polygon.
    return Polygon_type::CONCAVE;
  }

} // namespace internal
} // namespace Weights
} // namespace CGAL

#endif // CGAL_WEIGHTS_INTERNAL_POLYGON_UTILS_2_H
