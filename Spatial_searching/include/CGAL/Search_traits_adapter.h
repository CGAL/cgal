// Copyright (c) 2011 GeometryFactory (France).
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
// Author(s)     : Sebastien Loriot

#ifndef CGAL_SEARCH_TRAITS_WITH_INFO
#define CGAL_SEARCH_TRAITS_WITH_INFO

#include <CGAL/license/Spatial_searching.h>

#include <CGAL/disable_warnings.h>

#include <CGAL/Kd_tree_rectangle.h>
#include <CGAL/Euclidean_distance.h> //for default distance specialization
#include <CGAL/property_map.h>
#include <CGAL/assertions.h>

#include <boost/mpl/has_xxx.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/utility/enable_if.hpp>

namespace CGAL{

using ::get;
  
template <class Point_with_info,class PointPropertyMap,class Base_traits>
class Search_traits_adapter;
  
template <class Point_with_info,class PointPropertyMap,class Base_distance>
class Distance_adapter;
  
namespace internal{

BOOST_MPL_HAS_XXX_TRAIT_NAMED_DEF(Has_typedef_Iso_box_d,Iso_box_d,false)

template <class T,bool Has_iso_box_d=Has_typedef_Iso_box_d<T>::value >
struct Get_iso_box_d;

template <class T>
struct Get_iso_box_d<T,false>
{
  typedef void type;
};

template <class T>
struct Get_iso_box_d<T,true>
{
  typedef typename T::Iso_box_d type;
};
  
  template <class Point_with_info,class PointPropertyMap,class Base_traits>
  struct Spatial_searching_default_distance< ::CGAL::Search_traits_adapter<Point_with_info,PointPropertyMap,Base_traits> >{
    typedef ::CGAL::Distance_adapter<Point_with_info,
                                                 PointPropertyMap,
                                                 typename Spatial_searching_default_distance<Base_traits>::type> type;
  };

} //namespace internal
  
  
template <class Point_with_info,class PointPropertyMap,class Base_traits>
class Search_traits_adapter : public Base_traits{
  PointPropertyMap ppmap;

  CGAL_static_assertion( ( boost::is_same< boost::lvalue_property_map_tag,
                           typename boost::property_traits<PointPropertyMap>::category
                         >::value ) );
public:
  typedef Base_traits Base;
  typedef typename internal::Get_iso_box_d<Base>::type Iso_box_d;

  Search_traits_adapter(const PointPropertyMap& ppmap_=PointPropertyMap(),
                          const Base_traits& base=Base_traits()
  ):Base_traits(base),ppmap(ppmap_){}

  typedef typename Base_traits::Cartesian_const_iterator_d      Cartesian_const_iterator_d;
  typedef Point_with_info                                       Point_d;
  typedef typename Base_traits::FT                              FT;
  typedef typename Base_traits::Dimension                       Dimension;
  

  struct Construct_cartesian_const_iterator_d: public Base_traits::Construct_cartesian_const_iterator_d{
    PointPropertyMap ppmap;
    typedef typename Base_traits::Construct_cartesian_const_iterator_d Base;
    
    Construct_cartesian_const_iterator_d(const typename Base_traits::Construct_cartesian_const_iterator_d& base, const PointPropertyMap& ppmap_)
      :Base_traits::Construct_cartesian_const_iterator_d(base), ppmap(ppmap_){}
    
    typename Base_traits::Cartesian_const_iterator_d operator()(const Point_with_info& p) const
    { return Base::operator() (get(ppmap,p)); }

    typename Base_traits::Cartesian_const_iterator_d operator()(const Point_with_info& p, int)  const
    { return Base::operator() (get(ppmap,p),0); }

    // These 2 additional operators forward the call to Base_traits.
    // This is needed because of an undocumented requirement of 
    // Orthogonal_k_neighbor_search and Orthogonal_incremental_neighbor_search: 
    // Traits::Construct_cartesian_const_iterator should be callable 
    // on the query point type. If the query point type is the same as
    // Point_with_info, we disable it.

    template <typename Point> // boost::disable_if requires a template argument to work
    typename Base_traits::Cartesian_const_iterator_d operator()(const Point& p,
                                                                typename boost::disable_if<
                                                                boost::is_same<Point_with_info,
                                                                Point> >::type* = 0
      ) const
    { return Base::operator() (p); }

    template <typename Point> // boost::disable_if requires a template argument to work
    typename Base_traits::Cartesian_const_iterator_d operator()(const Point& p, int,
                                                                typename boost::disable_if<
                                                                boost::is_same<Point_with_info,
                                                                Point> >::type* = 0
      ) const
    { return Base::operator() (p,0); }
  };
  
  struct Construct_iso_box_d: public Base::Construct_iso_box_d{
    PointPropertyMap ppmap;
    typedef typename Base_traits::FT  FT; // needed for VC++, because otherwise it is taken from the private typedef of the base class
    typedef typename Base::Construct_iso_box_d Base_functor;

    Iso_box_d operator() () const {
      return Base_functor::operator() ();
    }
    Iso_box_d operator() (const Point_with_info& p, const Point_with_info& q) const
    {
      return Base_functor::operator() (get(ppmap,p),get(ppmap,q));
    }
  };
  
  const PointPropertyMap& point_property_map() const {return ppmap;}
  
  Construct_cartesian_const_iterator_d construct_cartesian_const_iterator_d_object() const {
    return Construct_cartesian_const_iterator_d(
      Base::construct_cartesian_const_iterator_d_object(),
      ppmap);
  }
};

template <class Point_with_info,class PointPropertyMap,class Base_distance>
class Distance_adapter : public Base_distance {
  PointPropertyMap ppmap;
  typedef typename Base_distance::FT FT;

  CGAL_static_assertion( ( boost::is_same< boost::lvalue_property_map_tag,
                           typename boost::property_traits<PointPropertyMap>::category
                         >::value ) );
public:
    
  Distance_adapter( const PointPropertyMap& ppmap_=PointPropertyMap(),
                         const Base_distance& distance=Base_distance()
  ):Base_distance(distance),ppmap(ppmap_){}

  using Base_distance::transformed_distance;
  
  typedef Point_with_info Point_d;
  typedef typename Base_distance::Query_item Query_item;

  const PointPropertyMap& point_property_map() const {return ppmap;}    

  FT transformed_distance(const Query_item& p1, const Point_with_info& p2) const
  {
    return Base_distance::transformed_distance(p1,get(ppmap,p2));
  }

  template <class FT,class Dimension>
  FT min_distance_to_rectangle(const Query_item& p, const CGAL::Kd_tree_rectangle<FT,Dimension>& b) const
  {
    return Base_distance::min_distance_to_rectangle(p,b);
  }

  template <class FT,class Dimension>
  FT min_distance_to_rectangle(const Query_item& p, const CGAL::Kd_tree_rectangle<FT,Dimension>& b,std::vector<FT>& dists) const
  {
    return Base_distance::min_distance_to_rectangle(p,b,dists);
  }

  template <class FT,class Dimension>
  FT max_distance_to_rectangle(const Query_item& p,const CGAL::Kd_tree_rectangle<FT,Dimension>& b) const
  {
    return Base_distance::max_distance_to_rectangle(p,b);
  }  
  template <class FT,class Dimension>
  FT max_distance_to_rectangle(const Query_item& p,const CGAL::Kd_tree_rectangle<FT,Dimension>& b,std::vector<FT>& dists) const
  {
    return Base_distance::max_distance_to_rectangle(p,b,dists);
  }  
};

}//namespace CGAL

#include <CGAL/enable_warnings.h>

#endif //CGAL_SEARCH_TRAITS_WITH_INFO
