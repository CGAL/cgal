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

#ifndef CGAL_SHAPE_DETECTION_REGION_GROWING_INTERNAL_PROPERTY_MAP_H
#define CGAL_SHAPE_DETECTION_REGION_GROWING_INTERNAL_PROPERTY_MAP_H

#include <CGAL/license/Shape_detection.h>

// STL includes.
#include <map>
#include <vector>
#include <memory>

// Boost includes.
#include <boost/iterator/transform_iterator.hpp>

// CGAL includes.
#include <CGAL/assertions.h>
#include <CGAL/property_map.h>
#include <CGAL/boost/graph/property_maps.h>

// Internal includes.
#include <CGAL/Shape_detection/Region_growing/internal/region_growing_traits.h>

namespace CGAL {
namespace Shape_detection {
namespace internal {

  template<
  typename ItemRange,
  typename PropertyMap>
  class Item_property_map {

  public:
    using key_type = std::size_t;
    using value_type = typename PropertyMap::value_type;
    using reference = const value_type&;
    using category = boost::lvalue_property_map_tag;

    Item_property_map(
      const ItemRange& item_range,
      const PropertyMap& property_map) :
    m_item_range(item_range),
    m_property_map(property_map)
    { }

    reference operator[](const key_type item_index) const {
      CGAL_precondition(item_index < m_item_range.size());
      const auto& key = *(m_item_range.begin() + item_index);
      return get(m_property_map, key);
    }

    friend inline reference get(const Item_property_map& item_map,
                                const key_type key) {
      return item_map[key];
    }

  private:
    const ItemRange& m_item_range;
    const PropertyMap& m_property_map;
  };

  template<typename Item>
  class Item_to_index_property_map {

  public:
    using key_type = Item;
    using value_type = std::size_t;
    using reference = value_type;
    using category = boost::readable_property_map_tag;

    template <class ItemRange>
    Item_to_index_property_map(const ItemRange& item_range)
    {
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

    friend inline value_type get(const Item_to_index_property_map& item_to_index_map,
                                 const key_type& key) {
      return item_to_index_map[key];
    }

  private:
    std::map<key_type, value_type> m_item_map;
  };

  class Seed_property_map {

  public:
    using key_type = std::size_t;
    using value_type = std::size_t;
    using category = boost::readable_property_map_tag;

    Seed_property_map(const std::vector<std::size_t>& seeds)
      : m_seeds(seeds)
    { }

    value_type operator[](const key_type key) const {
      return m_seeds[key];
    }

    friend value_type get(const Seed_property_map& seed_map,
                          const key_type key) {
      return seed_map[key];
    }

  private:
    const std::vector<std::size_t>& m_seeds;
  };

  template<typename ItemToIndexMap>
  class Item_to_region_index_map {

  public:
    using key_type = typename ItemToIndexMap::key_type;
    using value_type = long;
    using reference = value_type;
    using category = boost::readable_property_map_tag;

    Item_to_region_index_map() :
    m_item_to_index_map(nullptr)
    { }

    template<typename ItemRange, typename RegionRange>
    Item_to_region_index_map(
      const ItemRange& items,
      const ItemToIndexMap& item_to_index_map,
      const RegionRange& regions) :
    m_item_to_index_map(std::make_shared<ItemToIndexMap>(item_to_index_map)),
    m_indices(items.size(), -1) {

      long region_index = 0;
      for (const auto& region : regions) {
        for (const auto index : region.second) {
          CGAL_precondition(index < m_indices.size());
          m_indices[index] = region_index;
        }
        ++region_index;
      }
      CGAL_precondition(region_index == static_cast<long>(regions.size()));
    }

    inline friend value_type get(
      const Item_to_region_index_map& item_to_region_index_map,
      const key_type& key) {

      const auto& indices = item_to_region_index_map.m_indices;
      const auto& item_to_index_map = item_to_region_index_map.m_item_to_index_map;
      const std::size_t item_index = get(*item_to_index_map, key);
      CGAL_precondition(item_index < indices.size());
      return indices[item_index];
    }

  private:
    const std::shared_ptr<ItemToIndexMap> m_item_to_index_map;
    std::vector<long> m_indices;
  };

} // namespace internal
} // namespace Shape_detection
} // namespace CGAL

#endif // CGAL_SHAPE_DETECTION_REGION_GROWING_INTERNAL_PROPERTY_MAP_H
