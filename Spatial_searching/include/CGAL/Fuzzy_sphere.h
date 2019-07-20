// Copyright (c) 2002,2011 Utrecht University (The Netherlands).
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
// Author(s)     : Hans Tangelder (<hanst@cs.uu.nl>)


#ifndef CGAL_FUZZY_SPHERE_H
#define CGAL_FUZZY_SPHERE_H

#include <CGAL/license/Spatial_searching.h>

#include <CGAL/assertions.h>
#include <CGAL/Kd_tree_rectangle.h>
#include <CGAL/Search_traits_adapter.h>

#include <boost/type_traits.hpp>
#include <boost/utility/enable_if.hpp>

namespace CGAL {

namespace internal{

template <class SearchTraits,class Point_d>
class Fuzzy_sphere_impl
{
public:
  typedef typename SearchTraits::FT FT;
  typedef typename SearchTraits::Dimension Dimension;

private:
  SearchTraits traits;

  Point_d c;
  FT sq_radius;
  FT sq_inner_radius;
  FT sq_outer_radius;

public:
  // default constructor
  Fuzzy_sphere_impl(const SearchTraits& traits_=SearchTraits()):traits(traits_) {}

  // constructor
  Fuzzy_sphere_impl(const Point_d& center, FT r, FT eps=FT(0),
                    const SearchTraits& traits_ = SearchTraits())
    : traits(traits_), c(center), sq_radius(r*r)
  {
    CGAL_precondition(r >= 0);
    CGAL_precondition(eps >= 0);

    sq_inner_radius = (eps > r) ? FT(-1) : (r - eps) * (r - eps);
    sq_outer_radius = (r + eps) * (r + eps);
  }

  bool contains(const typename SearchTraits::Point_d& p) const {
    // test whether the distance between c and p is less than the radius
    FT distance=FT(0);
    typename SearchTraits::Construct_cartesian_const_iterator_d construct_it=
        traits.construct_cartesian_const_iterator_d_object();
    typename SearchTraits::Cartesian_const_iterator_d cit = construct_it(c),
        pit = construct_it(p), end = construct_it(c, 0);
    for (; cit != end && (distance <= sq_radius); ++cit, ++pit) {
      distance += CGAL::square((*cit)-(*pit));
    }

    return (distance <= sq_radius);
  }

  bool inner_range_intersects(const Kd_tree_rectangle<FT,Dimension>& rectangle) const {
    // test whether the sphere with radius (r-eps) intersects 'rectangle', i.e.
    // if the minimal distance of c to 'rectangle' is less than r-eps
    FT distance = FT(0);
    typename SearchTraits::Construct_cartesian_const_iterator_d construct_it=
        traits.construct_cartesian_const_iterator_d_object();
    typename SearchTraits::Cartesian_const_iterator_d cit = construct_it(c),
                                                      end = construct_it(c, 0);
    for (int i = 0; cit != end && (distance <= sq_inner_radius); ++cit, ++i) {
      if ((*cit) < rectangle.min_coord(i))
        distance += CGAL::square(rectangle.min_coord(i)-(*cit));
      else if ((*cit) > rectangle.max_coord(i))
        distance += CGAL::square((*cit)-rectangle.max_coord(i));
    }

    return (distance <= sq_inner_radius);
  }


  bool outer_range_contains(const Kd_tree_rectangle<FT,Dimension>& rectangle) const {
    // test whether the sphere with radius (r+eps) contains 'rectangle',
    // i.e. if the maximal distance of c to the boundary of 'rectangle'
    // is less than r+eps
    FT distance=FT(0);
    typename SearchTraits::Construct_cartesian_const_iterator_d construct_it=
        traits.construct_cartesian_const_iterator_d_object();
    typename SearchTraits::Cartesian_const_iterator_d cit = construct_it(c),
                                                      end = construct_it(c, 0);
    for (int i = 0; cit != end && (distance <= sq_outer_radius) ; ++cit,++i) {
      if ((*cit) <= (rectangle.min_coord(i)+rectangle.max_coord(i))/FT(2))
        distance += CGAL::square(rectangle.max_coord(i)-(*cit));
      else
        distance += CGAL::square((*cit)-rectangle.min_coord(i));
    }

    return (distance <= sq_outer_radius);
  }
}; // class Fuzzy_sphere_impl

}

template <class SearchTraits>
class Fuzzy_sphere:
    public internal::Fuzzy_sphere_impl<SearchTraits,typename SearchTraits::Point_d>
{
  typedef internal::Fuzzy_sphere_impl<SearchTraits,typename SearchTraits::Point_d> Base;
  typedef typename Base::FT FT;
public:
  // constructors
  Fuzzy_sphere(const SearchTraits& traits_=SearchTraits()):Base(traits_){}
  Fuzzy_sphere(const typename SearchTraits::Point_d& center, FT radius, FT epsilon=FT(0),const SearchTraits& traits_=SearchTraits()) :
    Base(center,radius,epsilon,traits_) {}
};

//specialization for Search_traits_adapter
template <class K,class PM,class Base_traits>
class Fuzzy_sphere< Search_traits_adapter<K,PM,Base_traits> >
  : public internal::Fuzzy_sphere_impl<Search_traits_adapter<K,PM,Base_traits>,typename Base_traits::Point_d>
{
  typedef Search_traits_adapter<K,PM,Base_traits> SearchTraits;
  typedef internal::Fuzzy_sphere_impl<SearchTraits,typename Base_traits::Point_d> Base;
  typedef typename Base_traits::FT FT;
public:
  // constructors
  Fuzzy_sphere(const SearchTraits& traits_=SearchTraits()):Base(traits_){}

  // Constructor for any point type that is not `SearchTraits::Point_d`
  template <typename Point> // boost::disable_if requires a template argument to work
  Fuzzy_sphere(const Point& center, FT radius, FT epsilon=FT(0),const SearchTraits& traits_=SearchTraits(),
               typename boost::disable_if<
               boost::is_same<typename SearchTraits::Point_d,
               Point> >::type* = 0)
    : Base(center,radius,epsilon,traits_) {}

  Fuzzy_sphere(const typename SearchTraits::Point_d& center, FT radius, FT epsilon=FT(0),
               const SearchTraits& traits_=SearchTraits())
    : Base(get(traits_.point_property_map(),center),radius,epsilon,traits_) {}
};

} // namespace CGAL
#endif // FUZZY_SPHERE_H
