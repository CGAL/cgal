// Copyright (c) 2011  GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Sebastien Loriot

#ifndef CGAL_SPATIAL_SORT_TRAITS_ADAPTER_2_H
#define CGAL_SPATIAL_SORT_TRAITS_ADAPTER_2_H

#include <CGAL/disable_warnings.h>

#include <boost/call_traits.hpp>

#include <CGAL/property_map.h>

namespace CGAL{

using ::get;

template<class Base_traits, class PointPropertyMap>
class Spatial_sort_traits_adapter_2
  : public Base_traits
{
  PointPropertyMap ppmap_;
public:
  Spatial_sort_traits_adapter_2(Base_traits base=Base_traits()):Base_traits(base){}

  Spatial_sort_traits_adapter_2(const PointPropertyMap& ppmap, Base_traits base=Base_traits())
    :Base_traits(base), ppmap_(ppmap){}

  typedef Base_traits Gt;
  typedef typename boost::property_traits<PointPropertyMap>::key_type Point_2;
  typedef typename boost::call_traits<Point_2>::param_type Arg_type;

  struct Less_x_2
    : public Base_traits::Less_x_2
  {
    Less_x_2(const PointPropertyMap& ppmap, const typename Base_traits::Less_x_2& base):
      Base_traits::Less_x_2(base), ppmap_(ppmap){}
    const PointPropertyMap& ppmap_;
    bool operator()(Arg_type p, Arg_type q) const {
      return static_cast<const typename Base_traits::Less_x_2*>(this)->operator()(get(ppmap_,p), get(ppmap_,q));
    }
  };

  struct Less_y_2
    : public Base_traits::Less_y_2
  {
    Less_y_2(const PointPropertyMap& ppmap, const typename Base_traits::Less_y_2& base):
      Base_traits::Less_y_2(base), ppmap_(ppmap){}
    const PointPropertyMap& ppmap_;
    bool operator()(Arg_type p, Arg_type q) const {
      return static_cast<const typename Base_traits::Less_y_2*>(this)->operator()(get(ppmap_,p), get(ppmap_,q));
    }
  };

  Less_x_2 less_x_2_object () const {return Less_x_2(ppmap_, static_cast<const Gt*>(this)->less_x_2_object() );}
  Less_y_2 less_y_2_object () const {return Less_y_2(ppmap_, static_cast<const Gt*>(this)->less_y_2_object() );}

  const PointPropertyMap& point_property_map() const {return ppmap_;}
};

} //namespace CGAL

#include <CGAL/enable_warnings.h>

#endif //CGAL_SPATIAL_SORT_TRAITS_ADAPTER_2_H
