// Copyright (c) 1998-2003  ETH Zurich (Switzerland).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
// 
//
// Author(s)     : Michael Hoffmann <hoffmann@inf.ethz.ch>

#ifndef CGAL_EXTREMAL_POLYGON_TRAITS_2_H
#define CGAL_EXTREMAL_POLYGON_TRAITS_2_H 1

#include <CGAL/license/Inscribed_areas.h>


#include <CGAL/Optimisation/assertions.h>
#include <CGAL/squared_distance_2.h>
#include <CGAL/Polygon_2.h>
#include <boost/bind.hpp>
#include <boost/function.hpp>

namespace CGAL {
template < class K_ >
struct Extremal_polygon_area_traits_2 {
  typedef          K_                                K;
  typedef typename K::FT                             FT;
  typedef typename K::Point_2                        Point_2;
#ifdef CGAL_OPTIMISATION_EXPENSIVE_PRECONDITION_TAG
  typedef typename K::Less_xy_2                      Less_xy_2;
  typedef typename K::Orientation_2                  Orientation_2;
#endif // CGAL_OPTIMISATION_EXPENSIVE_PRECONDITION_TAG

  struct Kgon_triangle_area {
    typedef FT              result_type;
  
    Kgon_triangle_area(const K& k_) : k(k_) {}
  
    result_type
    operator()(const Point_2& p, const Point_2& q, const Point_2& r) const
    { return CGAL_NTS abs(k.compute_area_2_object()(
      k.construct_triangle_2_object()(p, q, r))); }
  protected:
    K k;
  };

  typedef Kgon_triangle_area                         Baseop;
  typedef boost::function2<FT,Point_2,Point_2>       Operation;

  Extremal_polygon_area_traits_2() {}
  Extremal_polygon_area_traits_2(const K& k_) : k(k_) {}

  int min_k() const { return 3; }

  FT
  init( const Point_2&, const Point_2&) const
  { return FT( 0); }

  Operation
  operation( const Point_2& p) const
  { return boost::bind(Baseop(k), _1, _2, p); }

  template < class RandomAccessIC, class OutputIterator >
  OutputIterator
  compute_min_k_gon( RandomAccessIC points_begin,
                     RandomAccessIC points_end,
                     FT& max_area,
                     OutputIterator o) const
  // RandomAccessIC is a random access iterator or
  // circulator with value_type Point_2.
  // OutputIterator accepts int as value_type.
  //
  // PRE: n := | [points_begin, points_end) | >= min_k() and
  //  the points described by the range [points_begin, points_end)
  //  form the boundary of a convex polygon oriented counterclockwise.
  //
  // POST: write indices of the points from [points_begin, points_end)
  //  forming a min_k()-gon rooted at points_begin[0]
  //  of maximum area to o in counterclockwise order and return
  //  the past-the-end iterator for that range (== o + min_k()).
  {
    int number_of_points(
                         static_cast<int>(iterator_distance( points_begin, 
                                                             points_end)));
    CGAL_optimisation_precondition( number_of_points > min_k());
    
    // this gives the area of the triangle of two points with
    // the root:
    Operation op( operation( points_begin[0]));
    
    int p1( 1);
    int p2( 2);
    
    // maximal triangle so far (and the corresponding points):
    max_area = op( points_begin[p1], points_begin[p2]);
    int opt_p1( p1);
    int opt_p2( p2);
    
    // maximal triangle containing p1 so far:
    FT tmp_area( max_area);
    
    for (;;) {
      while ( p2 + 1 < number_of_points &&
              tmp_area < op( points_begin[p1], points_begin[p2+1])) {
        tmp_area = op( points_begin[p1], points_begin[++p2]);
      }
      if ( tmp_area > max_area) {
        max_area = tmp_area;
        opt_p1 = p1;
        opt_p2 = p2;
      }
      if ( ++p1 == number_of_points - 1)
        break;
      if ( p2 == p1)
        ++p2;
      tmp_area = op( points_begin[p1], points_begin[p2]);
    } // for (;;)
    
    // give result:
    *o++ = opt_p2;
    *o++ = opt_p1;
    *o++ = 0;
    return o;
  } // compute_min_k_gon( ... )

#ifdef CGAL_OPTIMISATION_EXPENSIVE_PRECONDITION_TAG
  Less_xy_2 less_xy_2_object() const { return k.less_xy_2_object(); }

  Orientation_2 orientation_2_object() const
  { return k.orientation_2_object(); }
#endif // CGAL_OPTIMISATION_EXPENSIVE_PRECONDITION_TAG

protected:
  K k;
};

} //namespace CGAL
#include <CGAL/Optimisation/assertions.h>
#include <cmath>
#ifdef CGAL_USE_LEDA
#include <CGAL/leda_real.h>
#endif // CGAL_USE_LEDA
namespace CGAL {

template < class K_ >
struct Extremal_polygon_perimeter_traits_2 {
  typedef          K_                                K;
  typedef typename K::FT                             FT;
  typedef typename K::Point_2                        Point_2;
#ifdef CGAL_OPTIMISATION_EXPENSIVE_PRECONDITION_TAG
  typedef typename K::Less_xy_2                    Less_xy_2;
  typedef typename K::Orientation_2                Orientation_2;
#endif // CGAL_OPTIMISATION_EXPENSIVE_PRECONDITION_TAG

  struct Kgon_triangle_perimeter {
    typedef FT              result_type;
  
    Kgon_triangle_perimeter(const K& k_): k(k_) {}

    // Added as workaround for VC2017 with /arch:AVX to fix
    // https://cgal.geometryfactory.com/CGAL/testsuite/CGAL-4.14-I-95/Inscribed_areas_Examples/TestReport_afabri_x64_Cygwin-Windows10_MSVC2017-Release-64bits.gz
    Kgon_triangle_perimeter(const Kgon_triangle_perimeter& other)
      : k(other.k)
    {}

    Kgon_triangle_perimeter& operator=(const Kgon_triangle_perimeter& other)
    {
      k = other.k;
      return *this;
    }
    
    result_type
    operator()(const Point_2& p, const Point_2& q, const Point_2& r) const
    { return dist(p, r) + dist(p, q) - dist(q, r); }
  
  private:
    result_type dist(const Point_2& p, const Point_2& q) const
    { return CGAL_NTS sqrt(k.compute_squared_distance_2_object()(p, q)); }
  
  protected:
    K k;
  };

  typedef Kgon_triangle_perimeter                    Baseop;
  typedef boost::function2<FT,Point_2,Point_2>       Operation;

  Extremal_polygon_perimeter_traits_2() {}
  Extremal_polygon_perimeter_traits_2(const K& k_) : k(k_) {}

  int min_k() const { return 2; }

  FT
  init( const Point_2& p, const Point_2& r) const
  { return operation( r)( p, r); }

  Operation
  operation( const Point_2& p) const
  { return boost::bind(Baseop(k), _1, _2, p); }

  template < class RandomAccessIC, class OutputIterator >
  OutputIterator
  compute_min_k_gon( RandomAccessIC points_begin,
                     RandomAccessIC points_end,
                     FT& max_perimeter,
                     OutputIterator o) const
  // RandomAccessIC is a random access iterator or
  // circulator with value_type Point_2.
  // OutputIterator accepts int as value_type.
  //
  // PRE: n := | [points_begin, points_end) | >= min_k() and
  //  the points described by the range [points_begin, points_end)
  //  form the boundary of a convex polygon oriented counterclockwise.
  //
  // POST: write indices of the points from [points_begin, points_end)
  //  forming a min_k()-gon rooted at points_begin[0] of maximum
  //  perimeter to o in counterclockwise order, set max_perimeter
  //  to twice this perimeter and return the past-the-end iterator
  //  for the range (== o + min_k()).
  {
    using std::less;
    using std::max_element;

    CGAL_optimisation_precondition_code(
      int number_of_points(
                           static_cast<int>(iterator_distance( points_begin, 
                                                               points_end)));)
    CGAL_optimisation_precondition( number_of_points > min_k());
    
    // kind of messy, but first we have to have something
    // like Distance (function object) ...
    RandomAccessIC maxi(
      max_element(
        points_begin + 1,
        points_end,
        boost::bind(
          less< FT >(),
          boost::bind(operation(points_begin[0]), _1, points_begin[0]),
          boost::bind(operation(points_begin[0]), _2, points_begin[0]))));
    
    // give result:
    max_perimeter = operation(*points_begin)(*maxi, *points_begin);
    *o++ = static_cast<int>(iterator_distance(points_begin, maxi));
    *o++ = 0;
    
    return o;
  } // compute_min_k_gon( ... )

#ifdef CGAL_OPTIMISATION_EXPENSIVE_PRECONDITION_TAG
  Less_xy_2 less_xy_2_object() const { return k.less_xy_2_object(); }

  Orientation_2 orientation_2_object() const
  { return k.orientation_2_object(); }
#endif // CGAL_OPTIMISATION_EXPENSIVE_PRECONDITION_TAG

protected:
  K k;
};


template < class RandomAccessIC, class OutputIterator >
inline
OutputIterator
maximum_area_inscribed_k_gon_2(
  RandomAccessIC points_begin,
  RandomAccessIC points_end,
  int k,
  OutputIterator o)
//
// preconditions:
// --------------
//  * Traits fulfills the requirements for an extremal polygon
//    traits class
//  * the range [points_begin, points_end) of size n > 0
//    describes the vertices of a convex polygon $P$
//    enumerated clock- or counterclockwise
//  * R is a CGAL representation class
//  * value_type of RandomAccessIC is Point_2<R>
//  * OutputIterator accepts Point_2<R> as value_type
//  * k >= 3
//
// functionality:
// --------------
// computes maximum area inscribed k-gon $P_k$
// of the polygon $P$,
// writes the indices (relative to points_begin)
// of $P_k$'s vertices to o and
// returns the past-the-end iterator of that sequence.
{
  typedef typename std::iterator_traits< RandomAccessIC >::value_type::R R;
  return extremal_polygon_2(
    points_begin,
    points_end,
    k,
    o,
    Extremal_polygon_area_traits_2< R >());
} // maximum_area_inscribed_k_gon_2( ... )

// backwards compatibility
template < class RandomAccessIC, class OutputIterator >
inline
OutputIterator
maximum_area_inscribed_k_gon(
  RandomAccessIC points_begin,
  RandomAccessIC points_end,
  int k,
  OutputIterator o)
{
  return maximum_area_inscribed_k_gon_2(
    points_begin, points_end, k, o);
} // maximum_area_inscribed_k_gon( ... )

template < class RandomAccessIC, class OutputIterator >
inline
OutputIterator
maximum_perimeter_inscribed_k_gon_2(
  RandomAccessIC points_begin,
  RandomAccessIC points_end,
  int k,
  OutputIterator o)
//
// preconditions:
// --------------
//  * Traits fulfills the requirements for an extremal polygon
//    traits class
//  * the range [points_begin, points_end) of size n > 0
//    describes the vertices of a convex polygon $P$
//    enumerated clock- or counterclockwise
//  * R is a CGAL representation class
//  * value_type of RandomAccessIC is Point_2<R>
//  * OutputIterator accepts Point_2<R> as value_type
//  * k >= 2
//
// functionality:
// --------------
// computes maximum perimeter inscribed k-gon $P_k$
// of the polygon $P$,
// writes the indices (relative to points_begin)
// of $P_k$'s vertices to o and
// returns the past-the-end iterator of that sequence.
{
  typedef typename std::iterator_traits< RandomAccessIC >::value_type::R R;
  return extremal_polygon_2(
    points_begin,
    points_end,
    k,
    o,
    Extremal_polygon_perimeter_traits_2< R >());
} // maximum_perimeter_inscribed_k_gon_2( ... )

// backwards compatibility
template < class RandomAccessIC, class OutputIterator >
inline
OutputIterator
maximum_perimeter_inscribed_k_gon(
  RandomAccessIC points_begin,
  RandomAccessIC points_end,
  int k,
  OutputIterator o)
{
  return maximum_perimeter_inscribed_k_gon_2(
    points_begin, points_end, k, o);
} // maximum_perimeter_inscribed_k_gon( ... )

} //namespace CGAL

#endif // ! (CGAL_EXTREMAL_POLYGON_TRAITS_2_H)
// ----------------------------------------------------------------------------
// ** EOF
// ----------------------------------------------------------------------------
