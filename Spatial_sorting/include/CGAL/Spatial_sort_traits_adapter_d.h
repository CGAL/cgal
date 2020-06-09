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


#ifndef CGAL_SPATIAL_SORT_TRAITS_ADAPTER_D_H
#define CGAL_SPATIAL_SORT_TRAITS_ADAPTER_D_H

#include <boost/call_traits.hpp>

#include <CGAL/property_map.h>

namespace CGAL{

using ::get;

template<class Base_traits, class PointPropertyMap>
class Spatial_sort_traits_adapter_d
  : public Base_traits
{
  PointPropertyMap ppmap_;
public:
  Spatial_sort_traits_adapter_d(Base_traits base=Base_traits()) : Base_traits(base){}

  Spatial_sort_traits_adapter_d(const PointPropertyMap& ppmap, Base_traits base=Base_traits())
    :Base_traits(base), ppmap_(ppmap){}

  typedef Base_traits Gt;
  typedef typename boost::property_traits<PointPropertyMap>::key_type Point_d;
  typedef typename boost::call_traits<Point_d>::param_type Arg_type;

  struct Point_dimension_d
    : public Base_traits::Point_dimension_d
  {
    Point_dimension_d(const PointPropertyMap& ppmap, const typename Base_traits::Point_dimension_d& base):
      Base_traits::Point_dimension_d(base), ppmap_(ppmap){}
    const PointPropertyMap& ppmap_;
    int operator()(Arg_type p) const {
      return static_cast<const typename Base_traits::Point_dimension_d*>(this)->operator()(get(ppmap_, p));
    }
  };

  struct Less_coordinate_d
    : public Base_traits::Less_coordinate_d
  {
    Less_coordinate_d(const PointPropertyMap& ppmap, const typename Base_traits::Less_coordinate_d& base):
      Base_traits::Less_coordinate_d(base), ppmap_(ppmap){}
    const PointPropertyMap& ppmap_;
    bool operator()(Arg_type p, Arg_type q, int i) const {
      return static_cast<const typename Base_traits::Less_coordinate_d*>(this)->operator()(get(ppmap_,p), get(ppmap_,q), i);
    }
  };

  struct Compute_coordinate_d
    : public Base_traits::Compute_coordinate_d
  {
    Compute_coordinate_d(const PointPropertyMap& ppmap, const typename Base_traits::Compute_coordinate_d& base):
      Base_traits::Compute_coordinate_d(base), ppmap_(ppmap){}
    const PointPropertyMap& ppmap_;
    bool operator()(Arg_type p, int i) const {
      return static_cast<const typename Base_traits::Compute_coordinate_d*>(this)->operator()(get(ppmap_,p), i);
    }
  };

  Point_dimension_d point_dimension_d_object () const {return Point_dimension_d(ppmap_, static_cast<const Gt*>(this)->point_dimension_d_object() );}
  Less_coordinate_d less_coordinate_d_object () const {return Less_coordinate_d(ppmap_, static_cast<const Gt*>(this)->less_coordinate_d_object() );}
  Compute_coordinate_d compute_coordinate_d_object () const {return Compute_coordinate_d(ppmap_, static_cast<const Gt*>(this)->compute_coordinate_d_object() );}

  const PointPropertyMap& point_property_map() const {return ppmap_;}
};

} //namespace CGAL

#endif //CGAL_SPATIAL_SORT_TRAITS_ADAPTER_D_H
