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

#ifndef CGAL_SHAPE_DETECTION_REGION_GROWING_SEGMENT_SET_SEGMENT_LENGTH_SORTING_H
#define CGAL_SHAPE_DETECTION_REGION_GROWING_SEGMENT_SET_SEGMENT_LENGTH_SORTING_H

#include <CGAL/license/Shape_detection.h>

#include <CGAL/Shape_detection/Region_growing/internal/property_map.h>
#include <CGAL/Shape_detection/Region_growing/internal/region_growing_traits.h>

namespace CGAL {
namespace Shape_detection {
namespace Segment_set {

/*!
  \ingroup PkgShapeDetectionRGOnSegments

  \brief Sorting of segments with respect to their length.

  `Items` of segments are sorted in decreasing length.

  \tparam GeomTraits
  a model of `Kernel`

  \tparam Item_
  a descriptor representing a given segment. Must be a model of `Hashable`.

  \tparam SegmentMap
  a model of `ReadablePropertyMap` whose key type is the value type of the input
  range and value type is `Kernel::Segment_2` or `Kernel::Segment_3`
*/
template<typename GeomTraits,
         typename Item_,
         typename SegmentMap>
class Segment_length_sorting
{
public:
  /// \name Types
  /// @{
  using Segment_map = SegmentMap;
  using Segment_type = typename boost::property_traits<Segment_map>::value_type;

  /// Item type.
  using Item = Item_;

  /// Seed range.
  using Seed_range = std::vector<Item>;
  /// @}

private:
  using FT = typename GeomTraits::FT;
  using Segment_set_traits = typename std::conditional<
    std::is_same<typename GeomTraits::Segment_2, Segment_type>::value,
    internal::Region_growing_traits_2<GeomTraits>,
    internal::Region_growing_traits_3<GeomTraits> >::type;

public:
  /// \name Initialization
  /// @{
  template<typename InputRange, typename NamedParameters = parameters::Default_named_parameters>
  Segment_length_sorting(const InputRange& input_range,
                         const NamedParameters& np = parameters::default_values())
  : m_segment_map(parameters::choose_parameter<SegmentMap>(parameters::get_parameter(
      np, internal_np::segment_map))),
      m_traits(parameters::choose_parameter<GeomTraits>(parameters::get_parameter(np, internal_np::geom_traits))),
      m_segment_set_traits(m_traits)
  {
    CGAL_precondition(input_range.size() > 0);

    using NP_helper = internal::Default_property_map_helper<NamedParameters, Item, typename InputRange::const_iterator, internal_np::item_map_t>;
    using Item_map = typename NP_helper::type;
    Item_map item_map = NP_helper::get(np);

    m_ordered.resize(input_range.size());
    std::size_t index = 0;
    for (auto it = input_range.begin(); it != input_range.end(); it++)
      m_ordered[index++] = get(item_map, it);
    m_scores.resize(input_range.size(), 0.);
  }
  /// @}

  /// \name Sorting
  /// @{
  void sort() {
    compute_scores();
    CGAL_precondition(m_scores.size() > 0);
    auto cmp = [this](const std::size_t i, const std::size_t j) {
      CGAL_precondition(i < m_scores.size());
      CGAL_precondition(j < m_scores.size());
      return m_scores[i] > m_scores[j];
    };
    std::vector<std::size_t> order(m_ordered.size());
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), cmp);
    std::vector<Item> tmp(m_ordered.size());
    for (std::size_t i = 0; i < m_ordered.size(); i++)
      tmp[i] = m_ordered[order[i]];
    m_ordered.swap(tmp);
  }
  /// @}

  /// \name Access
  /// @{
  const Seed_range &ordered() {
    return m_ordered;
  }
  /// @}

private:
  const Segment_map m_segment_map;
  const GeomTraits m_traits;
  const Segment_set_traits m_segment_set_traits;
  Seed_range m_ordered;
  std::vector<FT> m_scores;

  void compute_scores() {
    auto segment_length = m_segment_set_traits.compute_squared_length_object();
    std::size_t idx = 0;
    for (const Item& item : m_ordered) {
      const auto& segment = get(m_segment_map, item);
      FT len = segment_length(segment);
      m_scores[idx++] = CGAL::approximate_sqrt(len); // @todo no need to sqrt it?
    }
  }
};

} // namespace Segment_set
} // namespace Shape_detection
} // namespace CGAL

#endif // CGAL_SHAPE_DETECTION_REGION_GROWING_SEGMENT_SET_SEGMENT_LENGTH_SORTING_H
