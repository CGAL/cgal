// Copyright (c) 2018 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Florent Lafarge, Simon Giraudot, Thien Hoang, Dmitry Anisimov
//

#include <CGAL/license/Shape_detection.h>

// STL includes.
#include <map>
#include <vector>
#include <memory>

// CGAL includes.
#include <CGAL/assertions.h>
#include <CGAL/property_map.h>

#ifndef CGAL_SHAPE_DETECTION_REGION_GROWING_INTERNAL_PROPERTY_MAP_H
#define CGAL_SHAPE_DETECTION_REGION_GROWING_INTERNAL_PROPERTY_MAP_H

namespace CGAL {
namespace Shape_detection {
namespace internal {

  template<
  typename ItemRange,
  typename PropertyMap>
  class Item_property_map {

  public:
    using Item_range = ItemRange;
    using Property_map = PropertyMap;

    using value_type = typename Property_map::value_type;
    using reference = const value_type&;
    using key_type = std::size_t;
    using category = boost::lvalue_property_map_tag;

    Item_property_map(
      const Item_range& item_range,
      const Property_map& property_map) :
    m_item_range(item_range),
    m_property_map(property_map)
    { }

    reference operator[](const key_type item_index) const {
      CGAL_precondition(item_index < m_item_range.size());

      const auto& key = *(m_item_range.begin() + item_index);
      return get(m_property_map, key);
    }

    friend inline reference get(
      const Item_property_map& item_map,
      const key_type key) {

      return item_map[key];
    }

  private:
    const Item_range& m_item_range;
    const Property_map& m_property_map;
  };

  template<typename ItemRange>
  class Item_to_index_property_map {

  public:
    using Item_range = ItemRange;

    using Iterator = typename Item_range::const_iterator;
    using Item = typename Iterator::value_type;

    using value_type = std::size_t;
    using key_type = Item;
    using category = boost::lvalue_property_map_tag;

    using Item_map = std::map<key_type, value_type>;

    Item_to_index_property_map(const Item_range& item_range) :
    m_item_range(item_range) {

      value_type i = 0;
      for (const auto& item : item_range) {

        m_item_map[item] = i;
        ++i;
      }
    }

    value_type operator[](const key_type& key) const {

      const auto& value = m_item_map.find(key);

      if (value == m_item_map.end())
        return value_type(-1);

      return value->second;
    }

    friend inline value_type get(
      const Item_to_index_property_map& item_to_index_map,
      const key_type& key) {

      return item_to_index_map[key];
    }

  private:
    const Item_range& m_item_range;
    Item_map m_item_map;
  };

  class Seed_property_map {

  public:
    using key_type = std::size_t;
    using value_type = std::size_t;
    using category = boost::lvalue_property_map_tag;

    Seed_property_map(
      const std::vector<std::size_t>& seeds) :
    m_seeds(seeds)
    { }

    value_type operator[](const key_type key) const {
      return m_seeds[key];
    }

    friend value_type get(
      const Seed_property_map& seed_map,
      const key_type key) {

      return seed_map[key];
    }

  private:
    const std::vector<std::size_t>& m_seeds;
  };

} // namespace internal

namespace RG {

  class Point_to_shape_index_map {

  public:
    using key_type = std::size_t;
    using value_type = int;
    using reference = value_type;
    using category = boost::readable_property_map_tag;

    Point_to_shape_index_map() { }

    template<typename PointRange>
    Point_to_shape_index_map(
      const PointRange& points,
      const std::vector< std::vector<std::size_t> >& regions) :
    m_indices(new std::vector<int>(points.size(), -1)) {

      for (std::size_t i = 0; i < regions.size(); ++i)
        for (const std::size_t idx : regions[i])
          (*m_indices)[idx] = static_cast<int>(i);
    }

    inline friend value_type get(
      const Point_to_shape_index_map& point_to_shape_index_map,
      const key_type key) {

      const auto& indices = *(point_to_shape_index_map.m_indices);
      return indices[key];
    }

  private:
    std::shared_ptr< std::vector<int> > m_indices;
  };

} // namespace RG
} // namespace Shape_detection
} // namespace CGAL

#endif // CGAL_SHAPE_DETECTION_REGION_GROWING_INTERNAL_PROPERTY_MAP_H
