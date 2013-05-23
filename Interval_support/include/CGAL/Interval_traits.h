// Copyright (c) 2008 Max-Planck-Institute Saarbruecken (Germany).
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
//
// Author(s)     : Michael Hemmer <hemmer@mpi-inf.mpg.de>
//
// ============================================================================


/*! \file CGAL/Interval_traits.h
     This is experimental 
*/


/* bounds-related Interval functions */
// template<class Interval>  T lower(const Interval& x);
// template<class Interval>  T upper(const Interval& x);
// template<class Interval>  T width(const Interval& x);
// template<class Interval>  T median(const Interval& x);
// template<class Interval>  T norm(const Interval& x);

/* bounds-related Interval functions */
//// template<class Interval>  bool empty(const Interval& b); 
// template<class Interval>  bool singleton(const Interval& x);
// template<class Interval>  bool zero_in(const Interval& b);
// template<class Interval>  bool in(const T& r, const Interval& b);
// template<class Interval>  bool equal(const Interval& x, const Interval& y);
// template<class Interval>  bool overlap(const Interval& x, const Interval& y);
// template<class Interval>  bool subset(const Interval& a, const Interval& b);
// template<class Interval>  bool proper_subset(const Interval& a, const Interval& b);

/* set manipulation interval functions */
// template<class Interval>  Interval intersection(const Interval& x, const Interval& y);
// template<class Interval>  Interval hull(const Interval& x, const Interval& y);


#ifndef CGAL_INTERVAL_TRAITS_H
#define CGAL_INTERVAL_TRAITS_H

#include <CGAL/config.h>
#include <CGAL/tags.h>
#include <boost/type_traits/is_same.hpp>
#include <boost/utility/enable_if.hpp>

namespace CGAL {

namespace internal{

template<typename T> class Interval_traits_base{
public:
  typedef Interval_traits_base<T> Self; 
  typedef T                  Type;
  // typedef T                  Interval; 
  typedef CGAL::Tag_false    Is_interval; 
  typedef CGAL::Tag_false    With_empty_interval; 
  
  typedef CGAL::Null_functor Lower;
  typedef CGAL::Null_functor Upper; 
  typedef CGAL::Null_functor Width; 
  typedef CGAL::Null_functor Median;  
  typedef CGAL::Null_functor Norm; 
  typedef CGAL::Null_functor Empty;
  typedef CGAL::Null_functor Singleton;
  typedef CGAL::Null_functor In;
  typedef CGAL::Null_functor Zero_in;
  typedef CGAL::Null_functor Equal;
  typedef CGAL::Null_functor Overlap;
  typedef CGAL::Null_functor Subset;
  typedef CGAL::Null_functor Proper_Subset;
  typedef CGAL::Null_functor Intersection;
  typedef CGAL::Null_functor Hull;
};
}

template <typename T> class Interval_traits : public internal::Interval_traits_base<T>{};

class Exception_intersection_is_empty{}; 

// function returning type Bound 
template<typename Interval> inline 
typename Interval_traits<Interval>::Bound 
lower(const Interval& interval) {
    typename Interval_traits<Interval>::Lower lower;
    return lower(interval);
}

template<typename Interval> inline 
typename Interval_traits<Interval>::Bound 
upper(const Interval& interval) {
    typename Interval_traits<Interval>::Upper upper;
    return upper(interval);
}

template<typename Interval> inline 
typename Interval_traits<Interval>::Bound
width(Interval interval) {
    typename Interval_traits<Interval>::Width width;
    return width(interval);
}

template<typename Interval> inline 
typename Interval_traits<Interval>::Bound
median(Interval interval) {
    typename Interval_traits<Interval>::Median median;
    return median(interval);
}

template<typename Interval> inline 
typename Interval_traits<Interval>::Bound
norm(Interval interval) {
    typename Interval_traits<Interval>::Norm norm;
    return norm(interval);
}


// functions returning bool 

template<typename Interval> inline 
typename Interval_traits<Interval>::Empty::result_type 
empty(Interval interval) {
    typename Interval_traits<Interval>::Empty empty;
    return empty(interval);
}

template<typename Interval> inline 
typename Interval_traits<Interval>::Singleton::result_type  
singleton(Interval interval) {
    typename Interval_traits<Interval>::Singleton singleton;
    return singleton(interval);
}

template<typename Interval> inline 
typename Interval_traits<Interval>::In::result_type  
in(typename Interval_traits<Interval>::Bound x, Interval interval) {
    typename Interval_traits<Interval>::In in;
    return in(x,interval);
}

template<typename Interval> inline 
typename Interval_traits<Interval>::Zero_in::result_type
zero_in(Interval interval) {
    typename Interval_traits<Interval>::Zero_in zero_in;
    return zero_in(interval);
}

// This ones should be removed, since even boost_1_35_0 has changed to zero_in
template<typename Interval> inline 
typename Interval_traits<Interval>::Zero_in::result_type
in_zero(Interval interval) {
    typename Interval_traits<Interval>::Zero_in zero_in;
    return zero_in(interval);
}

template<typename Interval> inline 
typename Interval_traits<Interval>::Equal::result_type
equal(Interval interval1,Interval interval2) {
    typename Interval_traits<Interval>::Equal equal;
    return equal(interval1,interval2);
}

template<typename Interval> inline 
typename Interval_traits<Interval>::Overlap::result_type
overlap(Interval interval1, Interval interval2) {
    typename Interval_traits<Interval>::Overlap overlap;
    return overlap(interval1, interval2);
}

template<typename Interval> inline
typename Interval_traits<Interval>::Subset::result_type 
subset(Interval interval1, Interval interval2) {
    typename Interval_traits<Interval>::Subset subset;
    return subset(interval1, interval2);
}

template<typename Interval> inline
typename Interval_traits<Interval>::Proper_subset::result_type 
proper_subset(Interval interval1, Interval interval2) {
    typename Interval_traits<Interval>::Proper_subset proper_subset;
    return proper_subset(interval1, interval2);
}


// Set operations, functions returing Interval
//the enable_if is need for MSVC as it is not able to eliminate
//the function if Interval_traits<Interval>::Intersection has no result_type
//(like Null_functor)
template<typename Interval> inline 
typename Interval_traits<Interval>::Intersection::result_type
intersection(Interval interval1, Interval interval2, typename boost::enable_if<
             typename Interval_traits<Interval>::Is_interval
             >::type* = NULL
) {
    typename Interval_traits<Interval>::Intersection intersection;
    return intersection(interval1, interval2);
}

template<typename Interval> inline 
typename Interval_traits<Interval>::Hull::result_type
hull(Interval interval1, Interval interval2, typename boost::enable_if<
                                                       boost::is_same<
                                                          typename Interval_traits<Interval>::Is_interval,
                                                          Tag_true > >::type* = NULL) 
{
    typename Interval_traits<Interval>::Hull hull;
    return hull(interval1, interval2);
}



} //namespace CGAL

#endif // CGAL_INTERVAL_TRAITS_H
