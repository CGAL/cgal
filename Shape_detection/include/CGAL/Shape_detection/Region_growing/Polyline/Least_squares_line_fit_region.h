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

#ifndef CGAL_SHAPE_DETECTION_REGION_GROWING_POLYLINE_LEAST_SQUARES_LINE_FIT_REGION_H
#define CGAL_SHAPE_DETECTION_REGION_GROWING_POLYLINE_LEAST_SQUARES_LINE_FIT_REGION_H

#include <CGAL/license/Shape_detection.h>

// Boost includes.
#include <boost/unordered_map.hpp>

// Internal includes.
#include <CGAL/Shape_detection/Region_growing/internal/region_growing_traits.h>
#include <CGAL/use.h>

namespace CGAL {
namespace Shape_detection {
namespace Polyline {

  /*!
    \ingroup PkgShapeDetectionRGOnPolyline

    \brief Region type based on the quality of the least squares line
    fit applied to polyline vertices.

    This class fits a line, using \ref PkgPrincipalComponentAnalysisDRef "PCA",
    to chunks of polyline vertices and controls the quality of this fit.
    If all quality conditions are satisfied, the chunk is accepted as a valid region,
    otherwise rejected.

    \tparam GeomTraits
    a model of `Kernel`

    \tparam InputRange
    a model of `ConstRange` whose iterator type is `RandomAccessIterator`

    \tparam PointMap
    a model of `ReadablePropertyMap` whose key type is the value type of the input
    range and value type is `Kernel::Point_2` or `Kernel::Point_3`

    \cgalModels `RegionType`
  */
  template<
  typename GeomTraits,
  typename InputRange,
  typename PointMap>
  class Least_squares_line_fit_region {

  private:
    using Polyline_traits = typename std::conditional<
      std::is_same<typename GeomTraits::Point_2, typename PointMap::value_type>::value,
      internal::Region_growing_traits_2<GeomTraits>,
      internal::Region_growing_traits_3<GeomTraits> >::type;

  public:
    /// \name Types
    /// @{

    /// \cond SKIP_IN_MANUAL
    using Traits = GeomTraits;
    using Input_range = InputRange;
    using Point_map = PointMap;
    using Point_type = typename Point_map::value_type;
    /// \endcond

    /// Item type.
    using Item = typename InputRange::const_iterator;
    using Region = std::vector<Item>;

    /// Number type.
    typedef typename GeomTraits::FT FT;

    /// Primitive
    using Primitive = typename Polyline_traits::Line;
    using Result_type = std::vector<std::pair<Primitive, Region> >;

    /// Region map
    using Region_unordered_map = boost::unordered_map<Item, std::size_t, internal::hash_item<Item> >;
    using Region_index_map = boost::associative_property_map<Region_unordered_map>;
    /// @}

  private:
    using Point = typename Polyline_traits::Point;
    using Vector = typename Polyline_traits::Vector;
    using Line = typename Polyline_traits::Line;

    using Squared_length = typename Polyline_traits::Compute_squared_length;
    using Squared_distance = typename Polyline_traits::Compute_squared_distance;
    using Scalar_product = typename Polyline_traits::Compute_scalar_product;

  public:
    /// \name Initialization
    /// @{

    /*!
      \brief initializes all internal data structures.

      \tparam NamedParameters
      a sequence of \ref bgl_namedparameters "Named Parameters"

      \param input_range
      an instance of `InputRange` with polyline vertices

      \param np
      a sequence of \ref bgl_namedparameters "Named Parameters"
      among the ones listed below

      \cgalNamedParamsBegin
        \cgalParamNBegin{maximum_distance}
          \cgalParamDescription{the maximum distance from a vertex to a line}
          \cgalParamType{`GeomTraits::FT`}
          \cgalParamDefault{1}
        \cgalParamNEnd
        \cgalParamNBegin{maximum_angle}
          \cgalParamDescription{the maximum angle in degrees between
          the direction of a polyline edge and the direction of a line}
          \cgalParamType{`GeomTraits::FT`}
          \cgalParamDefault{25 degrees}
        \cgalParamNEnd
        \cgalParamNBegin{cosine_value}
          \cgalParamDescription{the cos value computed as `cos(maximum_angle * PI / 180)`,
          this parameter can be used instead of the `maximum_angle`}
          \cgalParamType{`GeomTraits::FT`}
          \cgalParamDefault{`cos(25 * PI / 180)`}
        \cgalParamNEnd
        \cgalParamNBegin{minimum_region_size}
          \cgalParamDescription{the minimum number of vertices a region must have}
          \cgalParamType{`std::size_t`}
          \cgalParamDefault{2}
        \cgalParamNEnd
        \cgalParamNBegin{point_map}
          \cgalParamDescription{an instance of `PointMap` that maps an item from `input_range`
          to `Kernel::Point_2` or `Kernel::Point_3`}
          \cgalParamDefault{`PointMap()`}
        \cgalParamNEnd
        \cgalParamNBegin{geom_traits}
          \cgalParamDescription{an instance of `GeomTraits`}
          \cgalParamDefault{`GeomTraits()`}
        \cgalParamNEnd
      \cgalNamedParamsEnd

      \pre `input_range.size() > 0`
      \pre `maximum_distance >= 0`
      \pre `maximum_angle >= 0 && maximum_angle <= 90`
      \pre `cosine_value >= 0 && cosine_value <= 1`
      \pre `minimum_region_size > 0`
    */
    template<typename NamedParameters = parameters::Default_named_parameters>
    Least_squares_line_fit_region(
      const InputRange& input_range,
      const NamedParameters& np = parameters::default_values()) :
    m_input_range(input_range),
    m_point_map(parameters::choose_parameter(parameters::get_parameter(
      np, internal_np::point_map), PointMap())),
    m_traits(parameters::choose_parameter(parameters::get_parameter(
      np, internal_np::geom_traits), GeomTraits())),
    m_polyline_traits(m_traits),
    m_squared_length(m_polyline_traits.compute_squared_length_object()),
    m_squared_distance(m_polyline_traits.compute_squared_distance_object()),
    m_scalar_product(m_polyline_traits.compute_scalar_product_object()) {

      CGAL_precondition(input_range.size() > 0);
      const FT max_distance = parameters::choose_parameter(
        parameters::get_parameter(np, internal_np::maximum_distance), FT(1));
      CGAL_precondition(max_distance >= FT(0));
      m_distance_threshold = max_distance;

      const FT max_angle = parameters::choose_parameter(
        parameters::get_parameter(np, internal_np::maximum_angle), FT(25));
      CGAL_precondition(max_angle >= FT(0) && max_angle <= FT(90));

      m_min_region_size = parameters::choose_parameter(
        parameters::get_parameter(np, internal_np::minimum_region_size), 2);
      CGAL_precondition(m_min_region_size > 0);

      const FT default_cos_value = static_cast<FT>(std::cos(CGAL::to_double(
        (max_angle * static_cast<FT>(CGAL_PI)) / FT(180))));
      const FT cos_value = parameters::choose_parameter(
        parameters::get_parameter(np, internal_np::cosine_value), default_cos_value);
      CGAL_precondition(cos_value >= FT(0) && cos_value <= FT(1));
      m_cos_value_threshold = cos_value;
    }

    /// @}

    /// \name Access
    /// @{

    /*!
      \brief implements `RegionType::region_index_map()`.

      This function creates an empty property map that maps iterators on the input range `Item` to std::size_t
    */
    Region_index_map region_index_map() {
      return Region_index_map(m_region_map);
    }

    /*!
      \brief implements `RegionType::primitive()`.

      This function provides the last primitive that has been fitted with the region.

      \return Primitive parameters that fits the region

      \pre `successful fitted primitive via successful call of update(region) with a sufficient large region`
    */
    Primitive primitive() const {
      return m_line_of_best_fit;
    }

    /*!
      \brief implements `RegionType::is_part_of_region()`.

      This function controls if a vertex with the item `query` is within
      the `maximum_distance` from the corresponding line and if the angle between the
      direction of the inward edge and the line's direction is within the `maximum_angle`.
      If both conditions are satisfied, it returns `true`, otherwise `false`.

      \param previous
      `Item` of the previous vertex

      \param query
      `Item` of the query vertex

      \return Boolean `true` or `false`

      \pre `region.size() > 0`
      \pre `previous` is a valid const_iterator of the input_range.
      \pre `query` is a valid const_iterator of the input_range.
    */
    bool is_part_of_region(
      const Item previous, const Item query,
      const Region&) {

      const Point& input_point = get(m_point_map, *previous);
      const Point& query_point = get(m_point_map, *query);

      // Update new reference line and direction.
      if (m_direction_of_best_fit == CGAL::NULL_VECTOR) {
        if (input_point == query_point) return true;
        CGAL_precondition(input_point != query_point);
        m_line_of_best_fit = Line(input_point, query_point);
        m_direction_of_best_fit = m_line_of_best_fit.to_vector();
        return true;
      }
      CGAL_precondition(m_direction_of_best_fit != CGAL::NULL_VECTOR);

      // Add equal points to the previously defined region.
      if (input_point == query_point) return true;
      CGAL_precondition(input_point != query_point);
      const Vector query_direction(input_point, query_point);

      // Check real conditions.
      const FT squared_distance_to_fitted_line =
        m_squared_distance(query_point, m_line_of_best_fit);
      const FT squared_distance_threshold =
        m_distance_threshold * m_distance_threshold;

      const FT cos_value =
        m_scalar_product(query_direction, m_direction_of_best_fit);
      const FT squared_cos_value = cos_value * cos_value;

      FT squared_cos_value_threshold =
        m_cos_value_threshold * m_cos_value_threshold;
      squared_cos_value_threshold *= m_squared_length(query_direction);
      squared_cos_value_threshold *= m_squared_length(m_direction_of_best_fit);

      return (
        ( squared_distance_to_fitted_line <= squared_distance_threshold ) &&
        ( squared_cos_value >= squared_cos_value_threshold ));
    }

    /*!
      \brief implements `RegionType::is_valid_region()`.

      This function controls if the `region` contains at least `minimum_region_size` vertices.

      \param region
      Vertices of the region represented as `Items`.

      \return Boolean `true` or `false`
    */
    inline bool is_valid_region(const Region& region) const {
      if (m_direction_of_best_fit == CGAL::NULL_VECTOR)
        return false; // all points are equal
      return (region.size() >= m_min_region_size);
    }

    /*!
      \brief implements `RegionType::update()`.

      This function fits the least squares line to all vertices from the `region`.

      \param region
      Vertices of the region represented as `Items`.

      \return Boolean `true` if the line fitting succeeded and `false` otherwise

      \pre `region.size() > 0`
    */
    bool update(const Region& region) {

      CGAL_precondition(region.size() > 0);
      if (region.size() == 1) { // create new reference line and direction
        m_direction_of_best_fit = CGAL::NULL_VECTOR;
      } else { // update reference line and direction
        if (m_direction_of_best_fit == CGAL::NULL_VECTOR)
          return false; // all points are equal
        CGAL_precondition(region.size() >= 2);
        std::tie(m_line_of_best_fit, m_direction_of_best_fit) =
          get_line_and_direction(region);
      }
      return true;
    }

    /// @}

    /// \cond SKIP_IN_MANUAL
    std::pair<Line, Vector> get_line_and_direction(
      const Region& region) const {

      // The best fit line will be a line fitted to all region points with
      // its direction being the line's direction.
      CGAL_precondition(region.size() > 0);
      const Line line_of_best_fit =
        m_polyline_traits.create_line(
          region, m_point_map).first;
      const Vector direction_of_best_fit =
        line_of_best_fit.to_vector();

      return std::make_pair(line_of_best_fit, direction_of_best_fit);
    }
    /// \endcond

  private:
    const Input_range& m_input_range;
    const Point_map m_point_map;
    const Traits m_traits;
    const Polyline_traits m_polyline_traits;
    Region_unordered_map m_region_map;

    FT m_distance_threshold;
    FT m_cos_value_threshold;
    std::size_t m_min_region_size;

    const Squared_length m_squared_length;
    const Squared_distance m_squared_distance;
    const Scalar_product m_scalar_product;

    Line m_line_of_best_fit;
    Vector m_direction_of_best_fit;
  };

} // namespace Polyline
} // namespace Shape_detection
} // namespace CGAL

#endif // CGAL_SHAPE_DETECTION_REGION_GROWING_POLYLINE_LEAST_SQUARES_LINE_FIT_REGION_H