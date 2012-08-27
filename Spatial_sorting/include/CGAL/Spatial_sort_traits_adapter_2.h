// Copyright (c) 2011  GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
//
// Author(s)     : Sebastien Loriot


#ifndef CGAL_SPATIAL_SORT_TRAITS_ADAPTER_2_H
#define CGAL_SPATIAL_SORT_TRAITS_ADAPTER_2_H

#include <boost/call_traits.hpp>
#include <boost/version.hpp>
#if BOOST_VERSION >= 104000
  #include <boost/property_map/property_map.hpp>
#else
  #include <boost/property_map.hpp>
#endif


namespace CGAL{
  
using ::get;

template<class Base_traits,class PointPropertyMap>
class Spatial_sort_traits_adapter_2:public Base_traits{
  PointPropertyMap ppmap_;
public:
  Spatial_sort_traits_adapter_2(Base_traits base=Base_traits()):Base_traits(base){}

  Spatial_sort_traits_adapter_2(const PointPropertyMap& ppmap,Base_traits base=Base_traits())
  :Base_traits(base),ppmap_(ppmap){}

  typedef Base_traits Gt;
  typedef typename boost::property_traits<PointPropertyMap>::key_type Point_2;
  typedef typename boost::call_traits<Point_2>::param_type Arg_type;

  struct Less_x_2 : public Base_traits::Less_x_2{
    Less_x_2(const PointPropertyMap& ppmap,const typename Base_traits::Less_x_2& base):
      Base_traits::Less_x_2(base),ppmap_(ppmap){}
    const PointPropertyMap& ppmap_;
    bool operator()(Arg_type p,Arg_type q) const {
      return static_cast<const typename Base_traits::Less_x_2*>(this)->operator()(get(ppmap_,p),get(ppmap_,q));
    }
  };

  struct Less_y_2 : public Base_traits::Less_y_2{
    Less_y_2(const PointPropertyMap& ppmap,const typename Base_traits::Less_y_2& base):
      Base_traits::Less_y_2(base),ppmap_(ppmap){}
    const PointPropertyMap& ppmap_;
    bool operator()(Arg_type p,Arg_type q) const {
      return static_cast<const typename Base_traits::Less_y_2*>(this)->operator()(get(ppmap_,p),get(ppmap_,q));
    }
  };

  Less_x_2 less_x_2_object () const {return Less_x_2(ppmap_,static_cast<const Gt*>(this)->less_x_2_object() );}
  Less_y_2 less_y_2_object () const {return Less_y_2(ppmap_,static_cast<const Gt*>(this)->less_y_2_object() );}
  
  const PointPropertyMap& point_property_map() const {return ppmap_;}

};

} //namespace CGAL

#endif //CGAL_SPATIAL_SORT_TRAITS_ADAPTER_2_H
