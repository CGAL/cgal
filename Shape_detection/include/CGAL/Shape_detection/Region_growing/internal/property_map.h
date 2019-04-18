// Copyright (c) 2018 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
//
//
// Author(s)     : Florent Lafarge, Simon Giraudot, Thien Hoang, Dmitry Anisimov
//

#include <CGAL/license/Shape_detection.h>

// STL includes.
#include <map>
#include <vector>

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
                
      CGAL_precondition(item_index >= 0);
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
} // namespace Shape_detection
} // namespace CGAL

#endif // CGAL_SHAPE_DETECTION_REGION_GROWING_INTERNAL_PROPERTY_MAP_H
