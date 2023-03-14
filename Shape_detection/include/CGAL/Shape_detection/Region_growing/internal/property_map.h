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

// CGAL includes.
#include <CGAL/property_map.h>
#include <CGAL/Named_function_parameters.h>
#include <CGAL/boost/graph/property_maps.h>
#include <CGAL/Default.h>

// Boost includes.
#include <boost/iterator/transform_iterator.hpp>

namespace CGAL {
namespace Shape_detection {
namespace internal {

template <class NP, class Item, class Iterator, class Tag,
          bool np_provided = /*true*/ !std::is_same<typename internal_np::Lookup_named_param_def<Tag,NP,Default>::type, Default>::value >
struct Default_property_map_helper
{
  using type = typename internal_np::Lookup_named_param_def<Tag,NP,Default>::type;
  static type get(const NP& np)
  {
    return parameters::get_parameter(np, Tag());
  }
};

template <class NP, class Item, class Iterator, class Tag>
struct Default_property_map_helper<NP, Item, Iterator, Tag, false>
{
  struct No_property_map_given_in_named_parameters_and_no_deducible_default {};
  static const bool iterator_is_item = std::is_same<Item, Iterator>::value;
  static const bool value_type_is_item = std::is_same<Item, typename std::iterator_traits<Iterator>::value_type>::value;
  using type = std::conditional_t<iterator_is_item,
                                  Identity_property_map<Item>,
                                  std::conditional_t<value_type_is_item,
                                                     Dereference_property_map<const Item, Iterator>,
                                                     No_property_map_given_in_named_parameters_and_no_deducible_default>>;

  static type get(const NP&)
  {
    return type();
  }
};

template <class ItemMap, class Item, class Iterator>
struct Item_map_helper
{
  using type = ItemMap;
  static const ItemMap& get(const ItemMap& m)
  {
    return m;
  }
};

template <class Item, class Iterator>
struct Item_map_helper<Default, Item, Iterator>
{
  using type = typename Default_property_map_helper<Default, Item, Iterator, int, false>::type;
  static type get(Default)
  {
    return type();
  }
};

} // namespace internal
} // namespace Shape_detection
} // namespace CGAL

#endif // CGAL_SHAPE_DETECTION_REGION_GROWING_INTERNAL_PROPERTY_MAP_H
