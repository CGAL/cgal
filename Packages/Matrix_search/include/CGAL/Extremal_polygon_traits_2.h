// ============================================================================
//
// Copyright (c) 1998 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       : $CGAL_Revision $
// release_date  : $CGAL_Date $
//
// file          : Extremal_polygon_traits_2.h
// chapter       : $CGAL_Chapter: Geometric Optimisation $
// package       : $CGAL_Package: Matrix_search $
// source        : mon_search.aw
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Michael Hoffmann <hoffmann@inf.ethz.ch>
//
// coordinator   : ETH Zurich (Bernd Gaertner <gaertner@inf.ethz.ch>)
//
// Predefined Traits classes for Extremal Polygon Computation
// ============================================================================

#if ! (CGAL_EXTREMAL_POLYGON_TRAITS_2_H)
#define CGAL_EXTREMAL_POLYGON_TRAITS_2_H 1

#include <CGAL/optimisation_assertions.h>
#ifndef CGAL_SQUARED_DISTANCE_2_H
#include <CGAL/squared_distance_2.h>
#endif // CGAL_SQUARED_DISTANCE_2_H
#ifndef CGAL_POLYGON_2_H
#include <CGAL/Polygon_2.h>
#endif // CGAL_POLYGON_2_H
#ifndef CGAL_FUNCTION_OBJECTS_H
#include <CGAL/function_objects.h>
#endif // CGAL_FUNCTION_OBJECTS_H

template < class R > inline
#ifndef CGAL_CFG_RETURN_TYPE_BUG_1
typename R::FT
#else
R_FT_return(R)
#endif
CGAL_Kgon_triangle_area( const CGAL_Point_2< R >& p,
                         const CGAL_Point_2< R >& q,
                         const CGAL_Point_2< R >& r)
{
  return CGAL_abs( p.x() * ( q.y() - r.y()) +
                   q.x() * ( r.y() - p.y()) +
                   r.x() * ( p.y() - q.y()));
}

template < class _R >
class CGAL__Kgon_area_operator
: public binary_function< CGAL_Point_2< _R >,
                          CGAL_Point_2< _R >,
                          typename _R::FT >
{
public:
  typedef _R                 R;
  typedef CGAL_Point_2< R >  Point_2;
  typedef typename R::FT     FT;

  CGAL__Kgon_area_operator( const Point_2& p)
  : root( p)
  {}

  FT
  operator()( const Point_2& p, const Point_2& q) const
  { return CGAL_Kgon_triangle_area( p, q, root); }

private:
  const Point_2& root;
};



template < class _R >
class CGAL_Kgon_area_traits
{
public:
  typedef          _R                   R;
  typedef          CGAL_Point_2< R >    Point_2;
  typedef typename _R::FT               FT;
  typedef CGAL__Kgon_area_operator< R > Operation;

  int
  min_k() const
  { return 3; }

  FT
  init( const Point_2&, const Point_2&) const
  { return FT( 0); }

  Operation
  operation( const Point_2& p) const
  { return Operation( p); }

#ifndef CGAL_CFG_NO_MEMBER_TEMPLATES
  template < class RandomAccessIC, class OutputIterator >
#else
  typedef typename vector< Point_2 >::iterator
    RandomAccessIC;
  typedef typename vector< int >::reverse_iterator
    OutputIterator;
#endif
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
  // POST: write the points of [points_begin, points_end)
  //  forming a min_k()-gon rooted at points_begin[0]
  //  of maximum area to o in counterclockwise order and return
  //  the past-the-end iterator for that range (== o + min_k()).
  {
    int number_of_points(
      CGAL_iterator_distance( points_begin, points_end));
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

  #ifndef CGAL_CFG_NO_MEMBER_TEMPLATES
  template < class RandomAccessIC >
  #endif
  bool
  is_convex( RandomAccessIC points_begin,
             RandomAccessIC points_end) const
  // PRE: value_type of RandomAccessIC is Point_2
  // POST: return true, iff the points [ points_begin, points_end)
  //   form a convex chain.
  {
    typedef CGAL_Polygon_traits_2< R >        P_traits;
    typedef vector< Point_2 >                 Cont;
    typedef CGAL_Polygon_2< P_traits, Cont >  Polygon_2;
  
    Polygon_2 p( points_begin, points_end);
    return p.is_convex();
  } // is_convex( points_begin, points_end)

};

#include <CGAL/optimisation_assertions.h>
#ifndef CGAL_PROTECT_MATH_H
#include <math.h>
#define CGAL_PROTECT_MATH_H
#endif // CGAL_PROTECT_MATH_H
#ifndef CGAL_LEDA_REAL_H
#include <CGAL/leda_real.h>
#endif // CGAL_LEDA_REAL_H

inline double
CGAL_sqrt( double x)
{ return sqrt( x); }

inline leda_real
CGAL_sqrt( leda_real x)
{ return sqrt( x); }

template < class _FT >
struct CGAL_Sqrt
: public binary_function< _FT, _FT, _FT >
{
  typedef _FT  FT;

  FT
  operator()( const FT& x) const
  { return CGAL_sqrt( x); }

};
template < class _R >
class CGAL__Kgon_perimeter_operator
: public binary_function< CGAL_Point_2< _R >,
                          CGAL_Point_2< _R >,
                          typename _R::FT >
{
public:
  typedef _R                 R;
  typedef CGAL_Point_2< R >  Point_2;
  typedef typename R::FT     FT;

  CGAL__Kgon_perimeter_operator( const Point_2& p)
  : root( p)
  {}

  FT
  operator()( const Point_2& p, const Point_2& q) const
  { return dist( p, root) + dist( p, q) - dist( q, root); }

private:
  static
  FT
  dist( const Point_2& p, const Point_2& q)
  { return CGAL_sqrt( CGAL_squared_distance( p, q)); }

  const Point_2& root;
};


template < class _R >
class CGAL_Kgon_perimeter_traits
{
public:
  typedef          _R                         R;
  typedef          CGAL_Point_2< R >          Point_2;
  typedef typename _R::FT                     FT;
  typedef CGAL__Kgon_perimeter_operator< R >  Operation;

  int
  min_k() const
  { return 2; }

  FT
  init( const Point_2& p, const Point_2& r) const
  { return operation( r)( p, r); }

  Operation
  operation( const Point_2& p) const
  { return Operation( p); }

#ifndef CGAL_CFG_NO_MEMBER_TEMPLATES
  template < class RandomAccessIC, class OutputIterator >
#else
  typedef typename vector< Point_2 >::iterator
    RandomAccessIC;
  typedef typename vector< int >::reverse_iterator
    OutputIterator;
#endif
  OutputIterator
  compute_min_k_gon( RandomAccessIC points_begin,
                     RandomAccessIC points_end,
                     FT& max_perimeter,
                     OutputIterator o) const
  // RandomAccessIC is a random access iterator or
  // circulator with value_type Point_2.
  // OutputIterator has value_type Point_2.
  //
  // PRE: n := | [points_begin, points_end) | >= min_k() and
  //  the points described by the range [points_begin, points_end)
  //  form the boundary of a convex polygon oriented counterclockwise.
  //
  // POST: write the points of [points_begin, points_end)
  //  forming a min_k()-gon rooted at points_begin[0] of maximum
  //  perimeter to o in counterclockwise order and return the
  //  past-the-end iterator for that range (== o + min_k()).
  {
    CGAL_optimisation_precondition_code(
      int number_of_points(
        CGAL_iterator_distance( points_begin, points_end));)
    CGAL_optimisation_precondition( number_of_points > min_k());
    
    // kind of messy, but first we have to have something
    // like CGAL_Distance (function object) ...
    RandomAccessIC maxi(
      max_element(
        points_begin + 1,
        points_end,
        CGAL_compose2_2(
          less< FT >(),
          bind2nd( operation( points_begin[0]), points_begin[0]),
          bind2nd( operation( points_begin[0]), points_begin[0]))));
    
    // give result:
    *o++ = CGAL_iterator_distance( points_begin, maxi);
    *o++ = 0;
    
    return o;
  } // compute_min_k_gon( ... )

  #ifndef CGAL_CFG_NO_MEMBER_TEMPLATES
  template < class RandomAccessIC >
  #endif
  bool
  is_convex( RandomAccessIC points_begin,
             RandomAccessIC points_end) const
  // PRE: value_type of RandomAccessIC is Point_2
  // POST: return true, iff the points [ points_begin, points_end)
  //   form a convex chain.
  {
    typedef CGAL_Polygon_traits_2< R >        P_traits;
    typedef vector< Point_2 >                 Cont;
    typedef CGAL_Polygon_2< P_traits, Cont >  Polygon_2;
  
    Polygon_2 p( points_begin, points_end);
    return p.is_convex();
  } // is_convex( points_begin, points_end)

};


template < class RandomAccessIC,
           class OutputIterator >
inline
OutputIterator
CGAL_maximum_area_inscribed_k_gon(
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
//  * value_type of RandomAccessIC (=: Point_2)
//    is CGAL_Point_2<R> for some representation class R
//  * OutputIterator accepts Point_2 as value_type
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
  return _CGAL_maximum_area_inscribed_k_gon(
    points_begin,
    points_end,
    k,
    o,
    value_type( points_begin));
} // CGAL_maximum_area_inscribed_k_gon( ... )

template < class RandomAccessIC,
           class OutputIterator,
           class R >
inline
OutputIterator
_CGAL_maximum_area_inscribed_k_gon(
  RandomAccessIC points_begin,
  RandomAccessIC points_end,
  int k,
  OutputIterator o,
  CGAL_Point_2< R >*)
//
// preconditions:
// --------------
//  * Traits fulfills the requirements for an extremal polygon
//    traits class
//  * the range [points_begin, points_end) of size n > 0
//    describes the vertices of a convex polygon $P$
//    enumerated clock- or counterclockwise
//  * R is a CGAL representation class
//  * value_type of RandomAccessIC is CGAL_Point_2<R>
//  * OutputIterator accepts CGAL_Point_2<R> as value_type
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
  return CGAL_extremal_polygon(
    points_begin,
    points_end,
    k,
    o,
    CGAL_Kgon_area_traits< R >());
} // _CGAL_maximum_area_inscribed_k_gon( ... )

template < class RandomAccessIC,
           class OutputIterator >
inline
OutputIterator
CGAL_maximum_perimeter_inscribed_k_gon(
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
//  * value_type of RandomAccessIC (=: Point_2)
//    is CGAL_Point_2<R> for some representation class R
//  * OutputIterator accepts Point_2 as value_type
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
  return _CGAL_maximum_perimeter_inscribed_k_gon(
    points_begin,
    points_end,
    k,
    o,
    value_type( points_begin));
} // CGAL_maximum_perimeter_inscribed_k_gon( ... )

template < class RandomAccessIC,
           class OutputIterator,
           class R >
inline
OutputIterator
_CGAL_maximum_perimeter_inscribed_k_gon(
  RandomAccessIC points_begin,
  RandomAccessIC points_end,
  int k,
  OutputIterator o,
  CGAL_Point_2< R >*)
//
// preconditions:
// --------------
//  * Traits fulfills the requirements for an extremal polygon
//    traits class
//  * the range [points_begin, points_end) of size n > 0
//    describes the vertices of a convex polygon $P$
//    enumerated clock- or counterclockwise
//  * R is a CGAL representation class
//  * value_type of RandomAccessIC is CGAL_Point_2<R>
//  * OutputIterator accepts CGAL_Point_2<R> as value_type
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
  return CGAL_extremal_polygon(
    points_begin,
    points_end,
    k,
    o,
    CGAL_Kgon_perimeter_traits< R >());
} // _CGAL_maximum_perimeter_inscribed_k_gon( ... )


#endif // ! (CGAL_EXTREMAL_POLYGON_TRAITS_2_H)

// ----------------------------------------------------------------------------
// ** EOF
// ----------------------------------------------------------------------------

