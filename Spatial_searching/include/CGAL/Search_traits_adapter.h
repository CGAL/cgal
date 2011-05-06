// Copyright (c) 2010 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
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
//
// Author(s)     : Sebastien Loriot

#ifndef CGAL_SEARCH_TRAITS_WITH_INFO
#define CGAL_SEARCH_TRAITS_WITH_INFO

#include <CGAL/Kd_tree_rectangle.h>
#include <boost/mpl/has_xxx.hpp>
#include <CGAL/Euclidean_distance.h> //for default distance specialization

namespace CGAL{

template <class Point_with_info,class Point_accessor,class Base_traits>
class Search_traits_adapter;
  
template <class Point_with_info,class Point_accessor,class Base_distance>
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
  
  template <class Point_with_info,class Point_accessor,class Base_traits>
  struct Spatial_searching_default_distance< ::CGAL::Search_traits_adapter<Point_with_info,Point_accessor,Base_traits> >{
    typedef ::CGAL::Distance_adapter<Point_with_info,
                                                 Point_accessor,
                                                 typename Spatial_searching_default_distance<Base_traits>::type> type;
  };

} //namespace internal
  
  
template <class Point_with_info,class Point_accessor,class Base_traits>
class Search_traits_adapter : public Base_traits{
  Point_accessor accessor;
public:
  typedef Base_traits Base;
  typedef typename internal::Get_iso_box_d<Base>::type Iso_box_d;

  Search_traits_adapter(const Point_accessor& accessor_=Point_accessor(),
                          const Base_traits& base=Base_traits()
  ):Base_traits(base),accessor(accessor_){}

  typedef typename Base_traits::Cartesian_const_iterator_d      Cartesian_const_iterator_d;
  typedef Point_with_info                                       Point_d;
  typedef typename Base_traits::FT                              FT;

  struct Construct_cartesian_const_iterator_d: public Base_traits::Construct_cartesian_const_iterator_d{
    Point_accessor accessor;
    using Base_traits::Construct_cartesian_const_iterator_d::operator();
    
    Construct_cartesian_const_iterator_d(const typename Base_traits::Construct_cartesian_const_iterator_d& base, const Point_accessor& accessor_)
      :Base_traits::Construct_cartesian_const_iterator_d(base), accessor(accessor_){}
    
    typename Base_traits::Cartesian_const_iterator_d operator()(const Point_with_info& p) const
    { return this->operator() (accessor[p]); }

    typename Base_traits::Cartesian_const_iterator_d operator()(const Point_with_info& p, int)  const
    { return this->operator() (accessor[p],0); }
  };
  
  struct Construct_iso_box_d: public Base::Construct_iso_box_d{
    Point_accessor accessor;
    
    Iso_box_d operator() () const {return this->Construct_iso_box_d();}
    Iso_box_d operator() (const Point_with_info& p, const Point_with_info& q, FT epsilon=FT(0)) const
    {
      return static_cast<const typename Base::Construct_iso_box_d* >(this)->operator() (accessor[p],accessor[q],epsilon);
    }
  };
  
  const Point_accessor& point_accessor() const {return accessor;}
  
  Construct_cartesian_const_iterator_d construct_cartesian_const_iterator_d_object() const {
    return Construct_cartesian_const_iterator_d(
      static_cast<const Base*>(this)->construct_cartesian_const_iterator_d_object(),
      accessor);
  }
};

template <class Point_with_info,class Point_accessor,class Base_distance>
class Distance_adapter : public Base_distance {
  Point_accessor accessor;
  typedef typename Base_distance::FT FT;
public:
    
  Distance_adapter( const Point_accessor& accessor_=Point_accessor(),
                         const Base_distance& distance=Base_distance()
  ):Base_distance(distance),accessor(accessor_){}

  using Base_distance::transformed_distance;
  
  typedef Point_with_info Point_d;
  typedef typename Base_distance::Query_item Query_item;

  const Point_accessor& point_accessor() const {return accessor;}    

  FT transformed_distance(const Query_item& p1, const Point_with_info& p2) const
  {
    return this->transformed_distance(p1,accessor[p2]);
  }

  template <class FT>
  FT min_distance_to_rectangle(const Query_item& p, const CGAL::Kd_tree_rectangle<FT>& b) const
  {
    return static_cast<const Base_distance*>(this)->min_distance_to_rectangle(p,b);
  }

  template <class FT>
  FT max_distance_to_rectangle(const Query_item& p,const CGAL::Kd_tree_rectangle<FT>& b) const
  {
    return static_cast<const Base_distance*>(this)->max_distance_to_rectangle(p,b);
  }  
};

}//namespace CGAL

#endif //CGAL_SEARCH_TRAITS_WITH_INFO
