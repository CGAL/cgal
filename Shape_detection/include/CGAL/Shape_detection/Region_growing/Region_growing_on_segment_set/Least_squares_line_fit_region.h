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

    \tparam InputRange
    a model of `ConstRange` whose iterator type is `RandomAccessIterator`

    \tparam SegmentMap
    a model of `LValuePropertyMap` whose key type is the value type of the input
    range and value type is `Kernel::Segment_2` or `Kernel::Segment_3`

    \cgalModels `RegionType`
  */
  template<
  typename GeomTraits,
  typename InputRange,
  typename SegmentMap>
  class Least_squares_line_fit_region {

  public:
    /// \name Types
    /// @{

    /// \cond SKIP_IN_MANUAL
    using Traits = GeomTraits;
    using Input_range = InputRange;
    using Segment_map = SegmentMap;
    using Segment_type = typename Segment_map::value_type;
    /// \endcond

    /// Number type.
    typedef typename GeomTraits::FT FT;

    /// @}

  private:
    using Segment_set_traits = typename std::conditional<
      std::is_same<typename Traits::Segment_2, Segment_type>::value,
      internal::Region_growing_traits_2<Traits>,
      internal::Region_growing_traits_3<Traits> >::type;

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

      \param input_range
      an instance of `InputRange` with 2D or 3D segments

      \param np
      a sequence of \ref bgl_namedparameters "Named Parameters"
      among the ones listed below

      \param segment_map
      an instance of `SegmentMap` that maps an item from `input_range`
      to `Kernel::Segment_2` or `Kernel::Segment_3`

      \param traits
      an instance of `GeomTraits`

      \cgalNamedParamsBegin
        \cgalParamNBegin{distance_threshold}
          \cgalParamDescription{the maximum distance from the furthest vertex of a segment to a line}
          \cgalParamType{`GeomTraits::FT`}
          \cgalParamDefault{1}
        \cgalParamNEnd
        \cgalParamNBegin{angle_threshold}
          \cgalParamDescription{the maximum accepted angle in degrees between
          the direction of a segment and the direction of a line}
          \cgalParamType{`GeomTraits::FT`}
          \cgalParamDefault{25 degrees}
        \cgalParamNEnd
        \cgalParamNBegin{cos_value_threshold}
          \cgalParamDescription{the cos value computed as `cos(angle_threshold * PI / 180)`,
          this parameter can be used instead of the `angle_threshold`}
          \cgalParamType{`GeomTraits::FT`}
          \cgalParamDefault{`cos(25 * PI / 180)`}
        \cgalParamNEnd
        \cgalParamNBegin{min_region_size}
          \cgalParamDescription{the minimum number of segments a region must have}
          \cgalParamType{`std::size_t`}
          \cgalParamDefault{1}
        \cgalParamNEnd
      \cgalNamedParamsEnd

      \pre `input_range.size() > 0`
      \pre `distance_threshold >= 0`
      \pre `angle_threshold >= 0 && angle_threshold <= 90`
      \pre `cos_value_threshold >= 0 && cos_value_threshold <= 1`
      \pre `min_region_size > 0`
    */
    template<typename NamedParameters>
    Least_squares_line_fit_region(
      const InputRange& input_range,
      const NamedParameters& np,
      const SegmentMap segment_map = SegmentMap(),
      const GeomTraits traits = GeomTraits()) :
    m_input_range(input_range),
    m_segment_map(segment_map),
    m_segment_set_traits(traits),
    m_squared_length(m_segment_set_traits.compute_squared_length_object()),
    m_squared_distance(m_segment_set_traits.compute_squared_distance_object()),
    m_scalar_product(m_segment_set_traits.compute_scalar_product_object()) {

      CGAL_precondition(input_range.size() > 0);
      m_distance_threshold = parameters::choose_parameter(
        parameters::get_parameter(np, internal_np::distance_threshold), FT(1));
      CGAL_precondition(m_distance_threshold >= FT(0));

      const FT angle_threshold = parameters::choose_parameter(
        parameters::get_parameter(np, internal_np::angle_threshold), FT(25));
      CGAL_precondition(angle_threshold >= FT(0) && angle_threshold <= FT(90));

      m_min_region_size = parameters::choose_parameter(
        parameters::get_parameter(np, internal_np::min_region_size), 1);
      CGAL_precondition(m_min_region_size > 0);

      const FT cos_value_threshold = static_cast<FT>(std::cos(CGAL::to_double(
        (angle_threshold * static_cast<FT>(CGAL_PI)) / FT(180))));
      m_cos_value_threshold = parameters::choose_parameter(
        parameters::get_parameter(np, internal_np::cos_value_threshold), cos_value_threshold);
      CGAL_precondition(m_cos_value_threshold >= FT(0) && m_cos_value_threshold <= FT(1));
    }

    /// @}

    /// \name Access
    /// @{

    /*!
      \brief implements `RegionType::is_part_of_region()`.

      This function controls if a segment with the index `query_index` is within
      the `distance_threshold` from the corresponding line and if the angle
      between the direction of this segment and the line's direction is within the `angle_threshold`.
      If both conditions are satisfied, it returns `true`, otherwise `false`.

      \param query_index
      index of the query segment

      The first and third parameters are not used in this implementation.

      \return Boolean `true` or `false`

      \pre `query_index < input_range.size()`
    */
    bool is_part_of_region(
      const std::size_t,
      const std::size_t query_index,
      const std::vector<std::size_t>&) const {

      CGAL_precondition(query_index < m_input_range.size());
      const auto& key = *(m_input_range.begin() + query_index);
      const Segment& query_segment = get(m_segment_map, key);
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

      This function controls if the `region` contains at least `min_region_size` segments.

      \param region
      indices of segments included in the region

      \return Boolean `true` or `false`
    */
    inline bool is_valid_region(const std::vector<std::size_t>& region) const {
      return (region.size() >= m_min_region_size);
    }

    /*!
      \brief implements `RegionType::update()`.

      This function fits the least squares line to all segments from the `region`.

      \param region
      indices of segments included in the region

      \return Boolean `true` if the line fitting succeeded and `false` otherwise

      \pre `region.size() > 0`
    */
    bool update(const std::vector<std::size_t>& region) {

      CGAL_precondition(region.size() > 0);
      if (region.size() == 1) { // create new reference line and direction
        const std::size_t segment_index = region[0];
        CGAL_precondition(segment_index < m_input_range.size());

        // The best fit line will be a line obtained from this segment
        // with the same direction.
        const auto& key = *(m_input_range.begin() + segment_index);
        const Segment& segment = get(m_segment_map, key);
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
      const std::vector<std::size_t>& region) const {

      // The best fit line will be a line fitted to all region segments with
      // its direction being the line's direction.
      CGAL_precondition(region.size() > 0);
      const Line line_of_best_fit =
        m_segment_set_traits.create_line(
          m_input_range, m_segment_map, region).first;
      const Vector direction_of_best_fit =
        line_of_best_fit.to_vector();

      return std::make_pair(line_of_best_fit, direction_of_best_fit);
    }
    /// \endcond

  private:
    const Input_range& m_input_range;

    FT m_distance_threshold;
    FT m_cos_value_threshold;
    std::size_t m_min_region_size;

    const Segment_map m_segment_map;
    const Segment_set_traits m_segment_set_traits;

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
