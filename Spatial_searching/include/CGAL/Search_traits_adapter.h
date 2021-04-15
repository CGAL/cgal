// Copyright (c) 2011 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
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
#include <boost/iterator/iterator_facade.hpp>

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

public:
  typedef Base_traits Base;
  typedef typename internal::Get_iso_box_d<Base>::type Iso_box_d;

  Search_traits_adapter(const PointPropertyMap& ppmap_=PointPropertyMap(),
                          const Base_traits& base=Base_traits()
  ):Base_traits(base),ppmap(ppmap_){}


  typedef Point_with_info                                       Point_d;
  typedef typename Base_traits::FT                              FT;
  typedef typename Base_traits::Dimension                       Dimension;


  // Default if point map is lvalue: use Construct_cartesian_const_iterator_d
  struct Construct_cartesian_const_iterator_d_lvalue: public Base_traits::Construct_cartesian_const_iterator_d{
    PointPropertyMap ppmap;
    typedef typename Base_traits::Construct_cartesian_const_iterator_d Base;

    Construct_cartesian_const_iterator_d_lvalue
    (const typename Base_traits::Construct_cartesian_const_iterator_d& base, const PointPropertyMap& ppmap_)
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

  // If point map is not lvalue, use this work-around that stores a
  // Point object in a shared pointer to avoid iterating on a temp
  // object
  class No_lvalue_iterator
    : public boost::iterator_facade<No_lvalue_iterator,
                                    typename std::iterator_traits<typename Base::Cartesian_const_iterator_d>::value_type,
                                    std::random_access_iterator_tag
                                     >
  {
    typedef boost::iterator_facade<No_lvalue_iterator,
                                   typename std::iterator_traits<typename Base::Cartesian_const_iterator_d>::value_type,
                                   std::random_access_iterator_tag
                                   > Facade;

    typedef typename std::iterator_traits<typename Base::Cartesian_const_iterator_d>::value_type
    Dereference_type;
    typedef typename boost::property_traits<PointPropertyMap>::value_type
    Point;

    boost::shared_ptr<Point> point;
    std::size_t idx;

  public:

    No_lvalue_iterator() : point(NULL), idx(0) { }
    No_lvalue_iterator(const Point& point) : point(new Point(point)), idx(0) { }
    No_lvalue_iterator(const Point& point, int) : point(new Point(point)), idx(Base::Dimension::value) { }

  private:

    friend class boost::iterator_core_access;
    void increment()
    {
      ++idx;
      CGAL_assertion(point != boost::shared_ptr<Point>());
    }
    void decrement()
    {
      --idx;
      CGAL_assertion(point != boost::shared_ptr<Point>());
    }

    void advance(std::ptrdiff_t n)
    {
      idx += n;
      CGAL_assertion(point != boost::shared_ptr<Point>());
    }

    std::ptrdiff_t distance_to(const No_lvalue_iterator& other) const
    {
      return other.idx - this->idx;

    }
    bool equal(const No_lvalue_iterator& other) const
    {
      return this->idx == other.idx;
    }

    Dereference_type&
    dereference() const
    {
      // Point::operator[] takes an int as parameter...
      return const_cast<Dereference_type&>((*point)[static_cast<int>(idx)]);
    }

  };

  // Alternative Construct_cartesian_const_iterator_d if the point map
  // is not lvalue (generates No_lvalue_iterator objects)
  struct Construct_cartesian_const_iterator_d_no_lvalue {
    typedef No_lvalue_iterator result_type;
    PointPropertyMap ppmap;

    Construct_cartesian_const_iterator_d_no_lvalue
    (const typename Base_traits::Construct_cartesian_const_iterator_d&, const PointPropertyMap& ppmap_)
      : ppmap(ppmap_) { }

    No_lvalue_iterator operator()(const Point_with_info& p) const
    { return No_lvalue_iterator(get(ppmap, p)); }

    No_lvalue_iterator operator()(const Point_with_info& p, int)  const
    { return No_lvalue_iterator(get(ppmap, p),0); }

    // These 2 additional operators forward the call to Base_traits.
    // This is needed because of an undocumented requirement of
    // Orthogonal_k_neighbor_search and Orthogonal_incremental_neighbor_search:
    // Traits::Construct_cartesian_const_iterator should be callable
    // on the query point type. If the query point type is the same as
    // Point_with_info, we disable it.

    template <typename Point> // boost::disable_if requires a template argument to work
    No_lvalue_iterator operator()(const Point& p,
                                  typename boost::disable_if<
                                  boost::is_same<Point_with_info,
                                  Point> >::type* = 0
      ) const
    { return No_lvalue_iterator(p); }

    template <typename Point> // boost::disable_if requires a template argument to work
    No_lvalue_iterator operator()(const Point& p, int,
                                  typename boost::disable_if<
                                  boost::is_same<Point_with_info,
                                  Point> >::type* = 0
      ) const
    { return No_lvalue_iterator(p,0); }
  };

  // Select type of iterator + construct class depending on whether
  // point map is lvalue or not
  typedef typename boost::mpl::if_
            <boost::is_same
              <boost::lvalue_property_map_tag,
               typename boost::property_traits<PointPropertyMap>::category >,
               typename Base::Cartesian_const_iterator_d,
               No_lvalue_iterator>::type
    Cartesian_const_iterator_d;
  typedef typename boost::mpl::if_
            <boost::is_same
              <boost::lvalue_property_map_tag,
               typename boost::property_traits<PointPropertyMap>::category >,
             Construct_cartesian_const_iterator_d_lvalue,
             Construct_cartesian_const_iterator_d_no_lvalue>::type
    Construct_cartesian_const_iterator_d;

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

public:

  Distance_adapter( const PointPropertyMap& ppmap_=PointPropertyMap(),
                         const Base_distance& distance=Base_distance()
  ):Base_distance(distance),ppmap(ppmap_){}

  using Base_distance::transformed_distance;

  typedef typename Base_distance::FT FT;
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
