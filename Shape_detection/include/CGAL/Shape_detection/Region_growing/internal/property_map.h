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
#include <CGAL/Named_function_parameters.h>

// Internal includes.
#include <CGAL/Shape_detection/Region_growing/internal/region_growing_traits.h>
#include <CGAL/Shape_detection/Region_growing/internal/utils.h>

namespace CGAL {
namespace Shape_detection {
namespace internal {

template<
  typename T,
  typename PropertyMap>
class Dereference_property_map_adaptor : public boost::put_get_helper<typename PropertyMap::reference, Dereference_property_map_adaptor< T, PropertyMap> > {
public:
  using key_type = T ;
  using value_type = typename PropertyMap::value_type;
  using reference = typename PropertyMap::reference;
  using category = boost::readable_property_map_tag;

  template <class U>
  reference operator[](U i) const {
    return get(m_property_map, *i);
  }
  Dereference_property_map_adaptor(const PropertyMap &property_map) : m_property_map(property_map) {}
  Dereference_property_map_adaptor() {}

  private:
    const PropertyMap m_property_map;
};

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
                                                     Dereference_property_map<const Item>,
                                                     No_property_map_given_in_named_parameters_and_no_deducible_default>>;

  static type get(const NP&)
  {
    return type();
  }
};



} // namespace internal
} // namespace Shape_detection
} // namespace CGAL

#endif // CGAL_SHAPE_DETECTION_REGION_GROWING_INTERNAL_PROPERTY_MAP_H
