// Copyright (c) 2014 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Dmitry Anisimov, David Bommes, Kai Hormann, Pierre Alliez
//

#ifndef CGAL_BARYCENTRIC_INTERNAL_UTILS_2_H
#define CGAL_BARYCENTRIC_INTERNAL_UTILS_2_H

#include <CGAL/license/Barycentric_coordinates_2.h>

// STL includes.
#include <set>
#include <map>
#include <list>
#include <tuple>
#include <sstream>
#include <fstream>
#include <cassert>

// CGAL headers.
#include <CGAL/array.h>
#include <CGAL/assertions.h>

// Boost headers.
#include <boost/mpl/has_xxx.hpp>
#include <optional>

// Internal includes.
#include <CGAL/Weights/internal/polygon_utils_2.h>
#include <CGAL/Barycentric_coordinates_2/barycentric_enum_2.h>

namespace CGAL {
namespace Barycentric_coordinates {
namespace internal {
  using namespace Weights::internal;

  enum class Query_point_location {

    // Query point is located at the vertex of the polygon.
    ON_VERTEX = 0,

    // Query point is located on the edge of the polygon.
    ON_EDGE = 1,

    // Query point is located in the polygon's interior.
    ON_BOUNDED_SIDE = 2,

    // Query point is located in the polygon's exterior.
    ON_UNBOUNDED_SIDE = 3,

    // Location is unspecified. Leads to all coordinates being set to zero.
    UNSPECIFIED = 4
  };

  template<typename GeomTraits>
  class Default_sqrt {

  private:
    using Traits = GeomTraits;
    using FT = typename Traits::FT;

  public:
    FT operator()(const FT value) const {
      return static_cast<FT>(
        CGAL::sqrt(CGAL::to_double(CGAL::abs(value))));
    }
  };

  BOOST_MPL_HAS_XXX_TRAIT_NAMED_DEF(Has_nested_type_Sqrt, Sqrt, false)

  // Case: do_not_use_default = false.
  template<typename GeomTraits,
  bool do_not_use_default = Has_nested_type_Sqrt<GeomTraits>::value>
  class Get_sqrt {

  public:
    using Traits = GeomTraits;
    using Sqrt = Default_sqrt<Traits>;

    static Sqrt sqrt_object(const Traits& ) {
      return Sqrt();
    }
  };

  // Case: do_not_use_default = true.
  template<typename GeomTraits>
  class Get_sqrt<GeomTraits, true> {

  public:
    using Traits = GeomTraits;
    using Sqrt = typename Traits::Sqrt;

    static Sqrt sqrt_object(const Traits& traits) {
      return traits.sqrt_object();
    }
  };

  // Get default values.
  template<typename OutputIterator>
  void get_default(
    const std::size_t n, OutputIterator output) {

    for (std::size_t i = 0; i < n; ++i) {
      *(output++) = 0;
    }
  }

  // Normalize values.
  template<typename FT>
  void normalize(std::vector<FT>& values) {

    FT sum = FT(0);
    for (const FT& value : values) {
      sum += value;
    }

    CGAL_assertion(sum != FT(0));
    if (sum == FT(0)) return;

    const FT inv_sum = FT(1) / sum;
    for (FT& value : values) {
      value *= inv_sum;
    }
  }

  // Compute barycentric coordinates along the line.
  template<
  typename OutputIterator,
  typename GeomTraits>
  OutputIterator linear_coordinates_2(
    const typename GeomTraits::Point_2& source,
    const typename GeomTraits::Point_2& target,
    const typename GeomTraits::Point_2& query,
    OutputIterator coordinates,
    const GeomTraits& traits) {

    CGAL_precondition(source != target);
    if (source == target) {
      get_default(2, coordinates);
      return coordinates;
    }

    // Number type.
    using FT = typename GeomTraits::FT;

    // Functions.
    const auto scalar_product_2   = traits.compute_scalar_product_2_object();
    const auto squared_distance_2 = traits.compute_squared_distance_2_object();

    // Project point onto segment.
    const FT opposite_scalar_product =
    scalar_product_2(query - target, source - target);

    // Compute coordinates.
    CGAL_assertion(source != target);
    const FT b0 = opposite_scalar_product / squared_distance_2(source, target);
    const FT b1 = FT(1) - b0;

    // Return coordinates.
    *(coordinates++) = b0;
    *(coordinates++) = b1;

    return coordinates;
  }

  // Compute barycentric coordinates in the plane.
  template<
  typename OutputIterator,
  typename GeomTraits>
  OutputIterator planar_coordinates_2(
    const typename GeomTraits::Point_2& p0,
    const typename GeomTraits::Point_2& p1,
    const typename GeomTraits::Point_2& p2,
    const typename GeomTraits::Point_2& query,
    OutputIterator coordinates,
    const GeomTraits& traits) {

    // Number type.
    using FT = typename GeomTraits::FT;

    // Functions.
    const auto area_2 = traits.compute_area_2_object();
    const FT total_area = area_2(p0, p1, p2);

    CGAL_precondition(total_area != FT(0));
    if (total_area == FT(0)) {
      get_default(3, coordinates);
      return coordinates;
    }

    // Compute some related sub-areas.
    const FT A1 = area_2(p1, p2, query);
    const FT A2 = area_2(p2, p0, query);

    // Compute the inverted total area of the triangle.
    CGAL_assertion(total_area != FT(0));
    const FT inverted_total_area = FT(1) / total_area;

    // Compute coordinates.
    const FT b0 = A1 * inverted_total_area;
    const FT b1 = A2 * inverted_total_area;
    const FT b2 = FT(1) - b0 - b1;

    // Return coordinates.
    *(coordinates++) = b0;
    *(coordinates++) = b1;
    *(coordinates++) = b2;

    return coordinates;
  }

  // Find a polygon edge that contains a query, if any (an approximate method).
  template<
  typename VertexRange,
  typename GeomTraits,
  typename PointMap>
  std::optional< std::pair<Query_point_location, std::size_t> >
  get_edge_index_approximate(
    const VertexRange& polygon,
    const typename GeomTraits::Point_2& query,
    const GeomTraits& traits,
    const PointMap point_map) {

    using FT = typename GeomTraits::FT;

    const auto cross_product_2 = traits.compute_determinant_2_object();
    const auto scalar_product_2 = traits.compute_scalar_product_2_object();
    const auto squared_distance_2 = traits.compute_squared_distance_2_object();
    const auto construct_vector_2 = traits.construct_vector_2_object();

    CGAL_precondition(polygon.size() >= 3);
    const std::size_t n = polygon.size();

    const FT half = FT(1) / FT(2);
    const FT tolerance = FT(1) / FT(100000);
    const FT sq_tolerance = tolerance * tolerance;

    for (std::size_t i = 0; i < n; ++i) {
      const auto& p1 = get(point_map, *(polygon.begin() + i));

      const FT sq_r = squared_distance_2(query, p1);
      if (sq_r < sq_tolerance) {
        return std::make_pair(Query_point_location::ON_VERTEX, i);
      }

      const std::size_t ip = (i + 1) % n;
      const auto& p2 = get(point_map, *(polygon.begin() + ip));

      const auto s1 = construct_vector_2(query, p1);
      const auto s2 = construct_vector_2(query, p2);

      const FT A = half * cross_product_2(s1, s2);
      const FT D = scalar_product_2(s1, s2);

      if (CGAL::abs(A) < tolerance && D < FT(0)) {
        return std::make_pair(Query_point_location::ON_EDGE, i);
      }
    }
    return std::nullopt;
  }

  // Why this one does not work for harmonic coordinates? - Due to the imprecisions in the Mesh_2 class.
  // Find a polygon edge that contains a query, if any (an exact method).
  template<
  typename VertexRange,
  typename GeomTraits,
  typename PointMap>
  std::optional< std::pair<Query_point_location, std::size_t> >
  get_edge_index_exact(
    const VertexRange& polygon,
    const typename GeomTraits::Point_2& query,
    const GeomTraits& traits,
    const PointMap point_map) {

    const auto collinear_2 = traits.collinear_2_object();
    const auto collinear_are_ordered_along_line_2 =
      traits.collinear_are_ordered_along_line_2_object();
    CGAL_precondition(polygon.size() >= 3);

    const std::size_t n = polygon.size();
    for (std::size_t i = 0; i < n; ++i) {
      const auto& p1 = get(point_map, *(polygon.begin() + i));

      if (p1 == query) {
        return std::make_pair(Query_point_location::ON_VERTEX, i);
      }

      const std::size_t ip = (i + 1) % n;
      const auto& p2 = get(point_map, *(polygon.begin() + ip));

      if (
        collinear_2(p1, p2, query) &&
        collinear_are_ordered_along_line_2(p1, query, p2)) {

        return std::make_pair(Query_point_location::ON_EDGE, i);
      }
    }
    return std::nullopt;
  }

  // Check whether a query point belongs to the last polygon edge.
  template<
  typename VertexRange,
  typename OutputIterator,
  typename GeomTraits,
  typename PointMap>
  std::pair<OutputIterator, bool> coordinates_on_last_edge_2(
    const VertexRange& polygon,
    const typename GeomTraits::Point_2& query,
    OutputIterator coordinates,
    const GeomTraits& traits,
    const PointMap point_map) {

    using FT = typename GeomTraits::FT;
    const std::size_t n = polygon.size();

    std::vector<FT> b;
    b.reserve(2);

    const std::size_t isource = n - 1;
    const std::size_t itarget = 0;

    const auto& source = get(point_map, *(polygon.begin() + isource));
    const auto& target = get(point_map, *(polygon.begin() + itarget));

    linear_coordinates_2(
      source, target, query, std::back_inserter(b), traits);
    *(coordinates++) = b[1];
    for (std::size_t i = 1; i < n - 1; ++i) {
      *(coordinates++) = FT(0);
    }
    *(coordinates++) = b[0];

    return std::make_pair(coordinates, true);
  }

  // Compute barycentric coordinates along the polygon boundary.
  template<
  typename VertexRange,
  typename OutputIterator,
  typename GeomTraits,
  typename PointMap>
  std::pair<OutputIterator, bool> boundary_coordinates_2(
    const VertexRange& polygon,
    const typename GeomTraits::Point_2& query,
    const Query_point_location location,
    const std::size_t index,
    OutputIterator coordinates,
    const GeomTraits& traits,
    const PointMap point_map) {

    using FT = typename GeomTraits::FT;
    const std::size_t n = polygon.size();

    // Compute coordinates with respect to the query point location.
    switch (location) {

      case Query_point_location::ON_VERTEX: {
        CGAL_assertion(index < n);

        for (std::size_t i = 0; i < n; ++i)
          if (i == index) {
            *(coordinates++) = FT(1);
          } else {
            *(coordinates++) = FT(0);
          }
        return std::make_pair(coordinates, true);
      }

      case Query_point_location::ON_EDGE: {
        CGAL_assertion(index < n);

        if (index == n - 1) {
          return coordinates_on_last_edge_2(
            polygon, query, coordinates, traits, point_map);
        }

        const std::size_t indexp = (index + 1) % n;
        const auto& source = get(point_map, *(polygon.begin() + index));
        const auto& target = get(point_map, *(polygon.begin() + indexp));

        for (std::size_t i = 0; i < n; ++i)
          if (i == index) {
            linear_coordinates_2(source, target, query, coordinates, traits); ++i;
          } else {
            *(coordinates++) = FT(0);
          }
        return std::make_pair(coordinates, true);
      }

      default: {
        internal::get_default(n, coordinates);
        return std::make_pair(coordinates, false);
      }
    }
    return std::make_pair(coordinates, false);
  }

  // Locate a point with respect to a polygon.
  template<
  typename VertexRange,
  typename GeomTraits,
  typename PointMap>
  std::optional< std::pair<Query_point_location, std::size_t> >
  locate_wrt_polygon_2(
    const VertexRange& polygon,
    const typename GeomTraits::Point_2& query,
    const GeomTraits& traits,
    const PointMap point_map) {

    const Edge_case type = bounded_side_2(
      polygon, query, traits, point_map);

    // Locate point with respect to different polygon locations.
    switch (type) {
      case Edge_case::INTERIOR:
        return std::make_pair(Query_point_location::ON_BOUNDED_SIDE, std::size_t(-1));
      case Edge_case::EXTERIOR:
        return std::make_pair(Query_point_location::ON_UNBOUNDED_SIDE, std::size_t(-1));
      case Edge_case::BOUNDARY:
        return get_edge_index_exact(polygon, query, traits, point_map);
      default:
        return std::make_pair(Query_point_location::UNSPECIFIED, std::size_t(-1));
    }
    return std::nullopt;
  }

} // namespace internal
} // namespace Barycentric_coordinates
} // namespace CGAL

#endif // CGAL_BARYCENTRIC_INTERNAL_UTILS_2_H
