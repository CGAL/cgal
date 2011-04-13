// ============================================================================
//
// Copyright (c) 1998, 1999, 2000 The CGAL Consortium
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
// maintainer    : Michael Hoffmann <hoffmann@inf.ethz.ch>
// coordinator   : ETH
//
// Predefined Traits classes for Extremal Polygon Computation
// ============================================================================

#if ! (CGAL_EXTREMAL_POLYGON_TRAITS_2_H)
#define CGAL_EXTREMAL_POLYGON_TRAITS_2_H 1

#include <CGAL/Optimisation/assertions.h>
#include <CGAL/squared_distance_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/function_objects.h>

CGAL_BEGIN_NAMESPACE
template < class R > inline
#ifndef CGAL_CFG_RETURN_TYPE_BUG_1
typename R::FT
#else
R_FT_return(R)
#endif
Kgon_triangle_area( const Point_2< R >& p,
                         const Point_2< R >& q,
                         const Point_2< R >& r)
{
  return CGAL_NTS abs( p.x() * ( q.y() - r.y()) +
                         q.x() * ( r.y() - p.y()) +
                         r.x() * ( p.y() - q.y()));
}

template < class R_ >
class Kgon_area_operator
: public CGAL_STD::binary_function< Point_2< R_ >,
                                    Point_2< R_ >,
                                    typename R_::FT >
{
public:
  typedef R_                 R;
  typedef Point_2< R >  Point_2;
  typedef typename R::FT     FT;

  Kgon_area_operator( const Point_2& p)
  : root( p)
  {}

  FT
  operator()( const Point_2& p, const Point_2& q) const
  { return Kgon_triangle_area( p, q, root); }

private:
  const Point_2& root;
};



template < class R_ >
class Extremal_polygon_area_traits_2
{
public:
  typedef          R_              R;
  typedef          Point_2< R >    Point_2;
  typedef typename R::FT           FT;
  typedef Kgon_area_operator< R >  Operation;

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
  typedef typename std::vector< Point_2 >::iterator
    RandomAccessIC;
  typedef typename std::vector< int >::reverse_iterator
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
  // POST: write indices of the points from [points_begin, points_end)
  //  forming a min_k()-gon rooted at points_begin[0]
  //  of maximum area to o in counterclockwise order and return
  //  the past-the-end iterator for that range (== o + min_k()).
  {
    int number_of_points(
      iterator_distance( points_begin, points_end));
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
    typedef Polygon_traits_2< R >        P_traits;
    typedef std::vector< Point_2 >       Cont;
    typedef Polygon_2< P_traits, Cont >  Polygon_2;
  
    Polygon_2 p( points_begin, points_end);
    return p.is_convex();
  } // is_convex( points_begin, points_end)

};

CGAL_END_NAMESPACE
#include <CGAL/Optimisation/assertions.h>
#include <cmath>
#ifdef CGAL_USE_LEDA
#include <CGAL/leda_real.h>
#endif // CGAL_USE_LEDA
CGAL_BEGIN_NAMESPACE

template < class FT_ >
struct Sqrt : public CGAL_STD::binary_function< FT_, FT_, FT_ >
{
  typedef FT_  FT;
  FT operator()(const FT& x) const { return CGAL_NTS sqrt(x); }
};
template < class R_ >
class Kgon_perimeter_operator
: public CGAL_STD::binary_function< Point_2< R_ >,
                                    Point_2< R_ >,
                                    typename R_::FT >
{
public:
  typedef R_              R;
  typedef Point_2< R >    Point_2;
  typedef typename R::FT  FT;

  Kgon_perimeter_operator( const Point_2& p)
  : root( p)
  {}

  FT
  operator()( const Point_2& p, const Point_2& q) const
  { return dist( p, root) + dist( p, q) - dist( q, root); }

private:
  static
  FT
  dist( const Point_2& p, const Point_2& q)
  { return CGAL_NTS sqrt( squared_distance( p, q)); }

  const Point_2& root;
};


template < class R_ >
class Extremal_polygon_perimeter_traits_2
{
public:
  typedef          R_                    R;
  typedef          Point_2< R >          Point_2;
  typedef typename R::FT                 FT;
  typedef Kgon_perimeter_operator< R >   Operation;

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
  typedef typename std::vector< Point_2 >::iterator
    RandomAccessIC;
  typedef typename std::vector< int >::reverse_iterator
    OutputIterator;
#endif
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
#ifndef CGAL_CFG_NO_NAMESPACE
    using std::bind2nd;
    using std::less;
    using std::max_element;
#endif

    CGAL_optimisation_precondition_code(
      int number_of_points(
        iterator_distance( points_begin, points_end));)
    CGAL_optimisation_precondition( number_of_points > min_k());
    
    // kind of messy, but first we have to have something
    // like Distance (function object) ...
    RandomAccessIC maxi(
      max_element(
        points_begin + 1,
        points_end,
        compose2_2(
          less< FT >(),
          bind2nd(operation(points_begin[0]), points_begin[0]),
          bind2nd(operation(points_begin[0]), points_begin[0]))));
    
    // give result:
    max_perimeter = operation(*points_begin)(*maxi, *points_begin);
    *o++ = iterator_distance(points_begin, maxi);
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
    typedef Polygon_traits_2< R >        P_traits;
    typedef std::vector< Point_2 >       Cont;
    typedef Polygon_2< P_traits, Cont >  Polygon_2;
  
    Polygon_2 p( points_begin, points_end);
    return p.is_convex();
  } // is_convex( points_begin, points_end)

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

CGAL_END_NAMESPACE

#endif // ! (CGAL_EXTREMAL_POLYGON_TRAITS_2_H)

// ----------------------------------------------------------------------------
// ** EOF
// ----------------------------------------------------------------------------

