// Copyright (c) 2025 GeometryFactory SARL (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Mael Rouxel-Labbé
//

#ifndef CGAL_SHAPE_DETECTION_REGION_GROWING_SEGMENT_SET_LINE_SEGMENT_REGION_H
#define CGAL_SHAPE_DETECTION_REGION_GROWING_SEGMENT_SET_LINE_SEGMENT_REGION_H

#include <CGAL/license/Shape_detection.h>

#include <CGAL/Shape_detection/Region_growing/internal/region_growing_traits.h>

namespace CGAL {
namespace Shape_detection {
namespace Segment_set {

/*!
  \ingroup PkgShapeDetectionRGOnSegments

  \brief Region type based on the line of the first segment selected.

  This class uses the supporting line of the first segment picked for the region
  and expands it for all segments with a direction close to that of the first segment
  (close being defined by the `maximum_distance` parameter) and such that endpoints are
  not far from that supporting line (far being defined by the `maximum_angle` or `cosine_of_maximum_angle` parameter).

  \tparam GeomTraits
  a model of `Kernel`

  \tparam Item_
  a descriptor representing a given segment. Must be a model of `Hashable`.

  \tparam SegmentMap
  a model of `ReadablePropertyMap` whose key type is `Item_`
  and value type is `GeomTraits::Segment_2` or `GeomTraits::Segment_3`

  \cgalModels{RegionType}
*/
template<typename GeomTraits,
         typename Item_,
         typename SegmentMap>
class Line_segment_region {
private:
  using Segment_set_traits = typename std::conditional<
    std::is_same<typename GeomTraits::Segment_2, typename SegmentMap::value_type>::value,
    internal::Region_growing_traits_2<GeomTraits>,
    internal::Region_growing_traits_3<GeomTraits> >::type;

public:
  /// \name Types
  /// @{
  using Segment_map = SegmentMap;
  using Segment_type = typename Segment_map::value_type;

  /// Number type.
  using FT = typename GeomTraits::FT;

  /// Item type.
  using Item = Item_;
  using Region = std::vector<Item>;
#ifdef DOXYGEN_RUNNING
  using Primitive = typename GeomTraits::Line_2 or typename GeomTraits::Line_3;
#else
  using Primitive = typename Segment_set_traits::Line;
#endif
  using Region_unordered_map = std::unordered_map<Item, std::size_t, internal::hash_item<Item>>;
  using Region_index_map = boost::associative_property_map<Region_unordered_map>;
  /// @}

private:
  using Point = typename Segment_set_traits::Point;
  using Vector = typename Segment_set_traits::Vector;
  using Line = typename Segment_set_traits::Line;
  using Squared_length = typename Segment_set_traits::Compute_squared_length;
  using Squared_distance = typename Segment_set_traits::Compute_squared_distance;
  using Scalar_product = typename Segment_set_traits::Compute_scalar_product;
  using Compare_squared_distance = typename Segment_set_traits::Compare_squared_distance;

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
        \cgalParamDescription{the cosine value `cos(maximum_angle * PI / 180)`, to be used instead of the parameter `maximum_angle()`}
        \cgalParamType{`GeomTraits::FT`}
        \cgalParamDefault{`cos(25 * PI / 180)`}
      \cgalParamNEnd
      \cgalParamNBegin{minimum_region_size}
        \cgalParamDescription{the minimum number of segments a region must have}
        \cgalParamType{`std::size_t`}
        \cgalParamDefault{1}
      \cgalParamNEnd
      \cgalParamNBegin{segment_map}
        \cgalParamDescription{an instance of `SegmentMap` that maps `Item_`
        to `GeomTraits::Segment_2` or `GeomTraits::Segment_3`}
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
  Line_segment_region(const NamedParameters& np = parameters::default_values())
    : m_segment_map(parameters::choose_parameter<SegmentMap>(parameters::get_parameter(
    np, internal_np::segment_map))),
      m_traits(parameters::choose_parameter<GeomTraits>(parameters::get_parameter(np, internal_np::geom_traits))),
      m_segment_set_traits(m_traits),
      m_squared_length(m_segment_set_traits.compute_squared_length_object()),
      m_squared_distance(m_segment_set_traits.compute_squared_distance_object()),
      m_scalar_product(m_segment_set_traits.compute_scalar_product_object()),
      m_compare_squared_distance(m_segment_set_traits.compare_squared_distance_object())
  {
    const FT max_distance = parameters::choose_parameter(
      parameters::get_parameter(np, internal_np::maximum_distance), FT(1));
    CGAL_precondition(max_distance >= FT(0));
    m_distance_threshold = max_distance;
    m_squared_distance_threshold = m_distance_threshold * m_distance_threshold;
    m_min_region_size = parameters::choose_parameter(
      parameters::get_parameter(np, internal_np::minimum_region_size), 1);
    CGAL_precondition(m_min_region_size > 0);
    const FT max_angle = parameters::choose_parameter(
      parameters::get_parameter(np, internal_np::maximum_angle), FT(25));
    const FT default_cos_value = static_cast<FT>(std::cos(CGAL::to_double(
      (max_angle * static_cast<FT>(CGAL_PI)) / FT(180))));
    const FT cos_value = parameters::choose_parameter(
      parameters::get_parameter(np, internal_np::cosine_of_maximum_angle), default_cos_value);
    CGAL_precondition(cos_value >= FT(0) && cos_value <= FT(1));
    m_cos_value_threshold = cos_value;
    m_squared_cos_value_threshold = m_cos_value_threshold * m_cos_value_threshold;
  }

  /// @}
  /// \name Access
  /// @{

  /*!
    \brief implements `RegionType::region_index_map()`.

    This function creates an empty property map that maps iterators on the input range `Item` to `std::size_t`.
  */
  Region_index_map region_index_map() {
    return Region_index_map(m_region_map);
  }

  /*!
    \brief implements `RegionType::primitive()`.

    This function provides the support line of the seed segment.

    \return Primitive parameters that fits the region

    \pre is_valid_region(region)
    \pre update(region)
  */
  Primitive primitive() const {
    // Compute the primitive (line) from the seed segment on the fly
    const auto& seed_seg = get(m_segment_map, m_seed);
    const Point& seed_s = seed_seg.source();
    const Point& seed_t = seed_seg.target();
    return Line(seed_s, seed_t);
  }

  /*!
    \brief implements `RegionType::is_part_of_region()`.

    This function determines if a segment with the index `query` is within
    the `maximum_distance` from the corresponding line and if the angle between the
    direction of this segment and the line's direction is within the `maximum_angle`.
    If both conditions are satisfied, it returns `true`, otherwise `false`.

    \param query item of the query segment
  */
  bool is_part_of_region(const Item query,
                         const Region&) const {
    const auto& seg = get(m_segment_map, query);
    const Point& s = seg.source();
    const Point& t = seg.target();
    const Vector dir(s, t);
    if (m_distance_threshold > FT(0)) {
      const auto& seed_seg = get(m_segment_map, m_seed);
      const Point& seed_s = seed_seg.source();
      const Point& seed_t = seed_seg.target();
      Line seed_line(seed_s, seed_t);
      if (m_compare_squared_distance(s, seed_line, m_squared_distance_threshold) != SMALLER ||
          m_compare_squared_distance(t, seed_line, m_squared_distance_threshold) != SMALLER)
        return false;
    }
    if (m_cos_value_threshold < FT(1)) {
      const auto& seed_seg = get(m_segment_map, m_seed);
      const Point& seed_s = seed_seg.source();
      const Point& seed_t = seed_seg.target();
      Vector seed_dir(seed_s, seed_t);
      FT cos_value = m_scalar_product(dir, seed_dir);
      FT squared_cos_value = cos_value * cos_value;
      FT threshold = m_squared_cos_value_threshold;
      threshold *= m_squared_length(dir);
      threshold *= m_squared_length(seed_dir);
      if (squared_cos_value < threshold)
        return false;
    }
    return true;
  }

  /*!
    \brief implements `RegionType::is_valid_region()`.

    This function controls if the `region` contains at least `minimum_region_size` segments.

    \param region
    Segments of the region represented as `Items`.
  */

  inline bool is_valid_region(const Region& region) const {
    return (region.size() >= m_min_region_size);
  }

  /*!
    \brief implements `RegionType::update()`.

    This function uses the support line of the seed segment as primitive.

    \param region
    Segments of the region represented as `Items`.

    \return Boolean `true` if the line fitting succeeded and `false` otherwise

    \pre `region.size() > 0`
  */

  bool update(const Region& region) {
    CGAL_precondition(region.size() > 0);
    if (region.size() == 1) {
      m_seed = region[0];
      const auto& seg = get(m_segment_map, m_seed);
      const Point& s = seg.source();
      const Point& t = seg.target();
      if (s == t) return false;
    }
    return true;
  }
  /// @}

private:
  const Segment_map m_segment_map;
  GeomTraits m_traits;
  Segment_set_traits m_segment_set_traits;
  Region_unordered_map m_region_map;
  FT m_distance_threshold;
  FT m_cos_value_threshold;
  std::size_t m_min_region_size;
  const Squared_length m_squared_length;
  const Squared_distance m_squared_distance;
  const Scalar_product m_scalar_product;
  const Compare_squared_distance m_compare_squared_distance;
  FT m_squared_distance_threshold;
  FT m_squared_cos_value_threshold;
  Item m_seed;
};

} // namespace Segment_set
} // namespace Shape_detection
} // namespace CGAL

#endif // CGAL_SHAPE_DETECTION_REGION_GROWING_SEGMENT_SET_LINE_SEGMENT_REGION_H
