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

#ifndef CGAL_SHAPE_DETECTION_REGION_GROWING_SEGMENT_SET_LEAST_SQUARES_LINE_FIT_REGION_H
#define CGAL_SHAPE_DETECTION_REGION_GROWING_SEGMENT_SET_LEAST_SQUARES_LINE_FIT_REGION_H

#include <CGAL/license/Shape_detection.h>

// CGAL includes.
#include <CGAL/Dynamic_property_map.h>

// Internal includes.
#include <CGAL/Shape_detection/Region_growing/internal/region_growing_traits.h>

namespace CGAL {
namespace Shape_detection {
namespace Segment_set {

  /*!
    \ingroup PkgShapeDetectionRGOnSegments

    \brief Region type based on the quality of the least squares line
    fit applied to a segment set.

    This class fits a line, using \ref PkgPrincipalComponentAnalysisDRef "PCA",
    to chunks of 2D or 3D segments and controls the quality of this fit.
    If all quality conditions are satisfied, the chunk is accepted as a valid region,
    otherwise rejected.

    \tparam GeomTraits
    a model of `Kernel`

    \tparam Item_
    a descriptor representing a given segment. Must be a model of `Hashable`.

    \tparam SegmentMap
    a model of `ReadablePropertyMap` whose key type is `Item`
    and value type is `Kernel::Segment_2` or `Kernel::Segment_3`

    \cgalModels{RegionType}
  */
  template<
  typename GeomTraits,
  typename Item_,
  typename SegmentMap>
  class Least_squares_line_fit_region {

  private:
    using Segment_set_traits = typename std::conditional<
      std::is_same<typename GeomTraits::Segment_2, typename SegmentMap::value_type>::value,
      internal::Region_growing_traits_2<GeomTraits>,
      internal::Region_growing_traits_3<GeomTraits> >::type;

  public:
    /// \name Types
    /// @{

    /// \cond SKIP_IN_MANUAL
    using Segment_map = SegmentMap;
    using Segment_type = typename Segment_map::value_type;
    /// \endcond

    /// Number type.
    typedef typename GeomTraits::FT FT;

    /// Item type.
    using Item = Item_;
    using Region = std::vector<Item>;

    /// Primitive type depends on the dimension of the input data.
#ifdef DOXYGEN_RUNNING
    using Primitive = typename GeomTraits::Line_2 or typename GeomTraits::Line_3
#else
    using Primitive = typename Segment_set_traits::Line;
#endif

    /// Region map
    using Region_unordered_map = std::unordered_map<Item, std::size_t, internal::hash_item<Item>>;
    using Region_index_map = boost::associative_property_map<Region_unordered_map>;
    /// @}

  private:
    using Point = typename Segment_set_traits::Point;
    using Segment = typename Segment_set_traits::Segment;
    using Vector = typename Segment_set_traits::Vector;
    using Line = typename Segment_set_traits::Line;

    using Squared_length = typename Segment_set_traits::Compute_squared_length;
    using Squared_distance = typename Segment_set_traits::Compute_squared_distance;
    using Scalar_product = typename Segment_set_traits::Compute_scalar_product;

  public:
    /// \name Initialization
    /// @{

    /*!
      \brief initializes all internal data structures.

      \tparam NamedParameters
      a sequence of \ref bgl_namedparameters "Named Parameters"

      \param np
      a sequence of \ref bgl_namedparameters "Named Parameters"
      among the ones listed below

      \tparam InputRange
      a model of `ConstRange` whose iterator type is `RandomAccessIterator`

      \cgalNamedParamsBegin
        \cgalParamNBegin{maximum_distance}
          \cgalParamDescription{the maximum distance from the furthest vertex of a segment to a line}
          \cgalParamType{`GeomTraits::FT`}
          \cgalParamDefault{1}
        \cgalParamNEnd
        \cgalParamNBegin{maximum_angle}
          \cgalParamDescription{the maximum angle in degrees between the direction of a segment and the direction of a line}
          \cgalParamType{`GeomTraits::FT`}
          \cgalParamDefault{25 degrees}
        \cgalParamNEnd
        \cgalParamNBegin{cosine_of_maximum_angle}
          \cgalParamDescription{the cosine value `cos(maximum_angle * PI / 180)` to be used instead of the parameter `maximum_angle()`}
          \cgalParamType{`GeomTraits::FT`}
          \cgalParamDefault{`cos(25 * PI / 180)`}
        \cgalParamNEnd
        \cgalParamNBegin{minimum_region_size}
          \cgalParamDescription{the minimum number of segments a region must have}
          \cgalParamType{`std::size_t`}
          \cgalParamDefault{1}
        \cgalParamNEnd
        \cgalParamNBegin{segment_map}
          \cgalParamDescription{an instance of `SegmentMap` that maps an item from `input_range`
          to `Kernel::Segment_2` or `Kernel::Segment_3`}
          \cgalParamDefault{`SegmentMap()`}
        \cgalParamNEnd
        \cgalParamNBegin{geom_traits}
          \cgalParamDescription{an instance of `GeomTraits`}
          \cgalParamDefault{`GeomTraits()`}
        \cgalParamNEnd
      \cgalNamedParamsEnd

      \pre `maximum_distance >= 0`
      \pre `maximum_angle >= 0 && maximum_angle <= 90`
      \pre `cosine_of_maximum_angle >= 0 && cosine_of_maximum_angle <= 1`
      \pre `minimum_region_size > 0`
    */
    template<typename NamedParameters = parameters::Default_named_parameters>
    Least_squares_line_fit_region(
      const NamedParameters& np = parameters::default_values()) :
    m_segment_map(parameters::choose_parameter<SegmentMap>(parameters::get_parameter(
      np, internal_np::segment_map))),
    m_traits(parameters::choose_parameter<GeomTraits>(parameters::get_parameter(np, internal_np::geom_traits))),
    m_segment_set_traits(m_traits),
    m_squared_length(m_segment_set_traits.compute_squared_length_object()),
    m_squared_distance(m_segment_set_traits.compute_squared_distance_object()),
    m_scalar_product(m_segment_set_traits.compute_scalar_product_object()) {

      const FT max_distance = parameters::choose_parameter(
        parameters::get_parameter(np, internal_np::maximum_distance), FT(1));
      CGAL_precondition(max_distance >= FT(0));
      m_distance_threshold = max_distance;

      const FT max_angle = parameters::choose_parameter(
        parameters::get_parameter(np, internal_np::maximum_angle), FT(25));
      CGAL_precondition(max_angle >= FT(0) && max_angle <= FT(90));

      m_min_region_size = parameters::choose_parameter(
        parameters::get_parameter(np, internal_np::minimum_region_size), 1);
      CGAL_precondition(m_min_region_size > 0);

      const FT default_cos_value = static_cast<FT>(std::cos(CGAL::to_double(
        (max_angle * static_cast<FT>(CGAL_PI)) / FT(180))));
      const FT cos_value = parameters::choose_parameter(
        parameters::get_parameter(np, internal_np::cosine_of_maximum_angle), default_cos_value);
      CGAL_precondition(cos_value >= FT(0) && cos_value <= FT(1));
      m_cos_value_threshold = cos_value;
    }

    /// @}

    /// \name Access
    /// @{

    /*!
      \brief implements `RegionType::region_index_map()`.

      This function creates an empty property map that maps iterators on the input range `Item` to std::size_t.
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

      This function controls if a segment with the index `query` is within
      the `maximum_distance` from the corresponding line and if the angle between the
      direction of this segment and the line's direction is within the `maximum_angle`.
      If both conditions are satisfied, it returns `true`, otherwise `false`.

      \param query
      `Item` of the query segment

      The last parameter is not used in this implementation.

      \return Boolean `true` or `false`

    */
    bool is_part_of_region(
      const Item query,
      const Region&) const {

      const Segment& query_segment = get(m_segment_map, query);
      const Point& query_source = query_segment.source();
      const Point& query_target = query_segment.target();
      const Vector query_direction(query_source, query_target);

      const FT squared_distance_to_fitted_line =
        get_max_squared_distance(query_segment);
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

      This function controls if the `region` contains at least `minimum_region_size` segments.

      \param region
      Segments of the region represented as `Items`.

      \return Boolean `true` or `false`
    */
    inline bool is_valid_region(const Region& region) const {
      return (region.size() >= m_min_region_size);
    }

    /*!
      \brief implements `RegionType::update()`.

      This function fits the least squares line to all segments from the `region`.

      \param region
      Segments of the region represented as `Items`.

      \return Boolean `true` if the line fitting succeeded and `false` otherwise

      \pre `region.size() > 0`
    */
    bool update(const Region& region) {

      CGAL_precondition(region.size() > 0);
      if (region.size() == 1) { // create new reference line and direction
        const Item item = region[0];

        // The best fit line will be a line obtained from this segment
        // with the same direction.
        const Segment& segment = get(m_segment_map, item);
        const Point& source = segment.source();
        const Point& target = segment.target();
        if (source == target) return false;

        CGAL_precondition(source != target);
        m_line_of_best_fit = Line(source, target);
        m_direction_of_best_fit = m_line_of_best_fit.to_vector();

      } else { // update reference line and direction
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

      // The best fit line will be a line fitted to all region segments with
      // its direction being the line's direction.
      CGAL_precondition(region.size() > 0);
      const Line line_of_best_fit =
        m_segment_set_traits.create_line(
          region, m_segment_map).first;
      const Vector direction_of_best_fit =
        line_of_best_fit.to_vector();

      return std::make_pair(line_of_best_fit, direction_of_best_fit);
    }
    /// \endcond

  private:
    const Segment_map m_segment_map;
    const GeomTraits m_traits;
    const Segment_set_traits m_segment_set_traits;
    Region_unordered_map m_region_map;

    FT m_distance_threshold;
    FT m_cos_value_threshold;
    std::size_t m_min_region_size;

    const Squared_length m_squared_length;
    const Squared_distance m_squared_distance;
    const Scalar_product m_scalar_product;

    Line m_line_of_best_fit;
    Vector m_direction_of_best_fit;

    // The maximum squared distance from the vertices of the segment
    // to the best fit line.
    FT get_max_squared_distance(const Segment& segment) const {

      const Point& source = segment.source();
      const Point& target = segment.target();
      const FT squared_distance_source =
        m_squared_distance(source, m_line_of_best_fit);
      const FT squared_distance_target =
        m_squared_distance(target, m_line_of_best_fit);
      return (CGAL::max)(squared_distance_source, squared_distance_target);
    }
  };

} // namespace Segment_set
} // namespace Shape_detection
} // namespace CGAL

#endif // CGAL_SHAPE_DETECTION_REGION_GROWING_SEGMENT_SET_LEAST_SQUARES_LINE_FIT_REGION_H
