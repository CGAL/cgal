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

} // namespace internal
} // namespace Shape_detection
} // namespace CGAL

#endif // CGAL_SHAPE_DETECTION_REGION_GROWING_INTERNAL_PROPERTY_MAP_H
