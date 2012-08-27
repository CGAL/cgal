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


#ifndef CGAL_SPATIAL_SORT_TRAITS_ADAPTER_D_H
#define CGAL_SPATIAL_SORT_TRAITS_ADAPTER_D_H

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
class Spatial_sort_traits_adapter_d:public Base_traits{
  PointPropertyMap ppmap_;
public:
  Spatial_sort_traits_adapter_d(Base_traits base=Base_traits()):Base_traits(base){}

  Spatial_sort_traits_adapter_d(const PointPropertyMap& ppmap,Base_traits base=Base_traits())
  :Base_traits(base),ppmap_(ppmap){}

  typedef Base_traits Gt;
  typedef typename boost::property_traits<PointPropertyMap>::key_type Point_d;
  typedef typename boost::call_traits<Point_d>::param_type Arg_type;

    
    
  struct Point_dimension_d: public Base_traits::Point_dimension_d{
    Point_dimension_d(const PointPropertyMap& ppmap,const typename Base_traits::Point_dimension_d& base):
      Base_traits::Point_dimension_d(base),ppmap_(ppmap){}
    const PointPropertyMap& ppmap_;
    int operator()(Arg_type p) const {
      return static_cast<const typename Base_traits::Point_dimension_d*>(this)->operator()(get(ppmap_,p));
    }
  };

  struct Less_coordinate_d: public Base_traits::Less_coordinate_d{
    Less_coordinate_d(const PointPropertyMap& ppmap,const typename Base_traits::Less_coordinate_d& base):
      Base_traits::Less_coordinate_d(base),ppmap_(ppmap){}
    const PointPropertyMap& ppmap_;
    bool operator()(Arg_type p,Arg_type q,int i) const {
      return static_cast<const typename Base_traits::Less_coordinate_d*>(this)->operator()(get(ppmap_,p),get(ppmap_,q),i);
    }
  };


  struct Compute_coordinate_d: public Base_traits::Compute_coordinate_d{
    Compute_coordinate_d(const PointPropertyMap& ppmap,const typename Base_traits::Compute_coordinate_d& base):
      Base_traits::Compute_coordinate_d(base),ppmap_(ppmap){}
    const PointPropertyMap& ppmap_;
    bool operator()(Arg_type p,int i) const {
      return static_cast<const typename Base_traits::Compute_coordinate_d*>(this)->operator()(get(ppmap_,p),i);
    }
  };
 
 

  Point_dimension_d point_dimension_d_object () const {return Point_dimension_d(ppmap_,static_cast<const Gt*>(this)->point_dimension_d_object() );}
  Less_coordinate_d less_coordinate_d_object () const {return Less_coordinate_d(ppmap_,static_cast<const Gt*>(this)->less_coordinate_d_object() );}
  Compute_coordinate_d compute_coordinate_d_object () const {return Compute_coordinate_d(ppmap_,static_cast<const Gt*>(this)->compute_coordinate_d_object() );}
  
  const PointPropertyMap& point_property_map() const {return ppmap_;}
  
};

} //namespace CGAL

#endif //CGAL_SPATIAL_SORT_TRAITS_ADAPTER_D_H
