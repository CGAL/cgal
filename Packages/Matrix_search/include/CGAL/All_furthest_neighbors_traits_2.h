#line 1404 "mon_search.aw"
#line 18 "code_formatting.awi"
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
// file          : All_furthest_neighbors_traits_2.h
// chapter       : $CGAL_Chapter: Geometric Optimisation $
// package       : $CGAL_Package: Matrix_search $
// source        : mon_search.aw
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Michael Hoffmann <hoffmann@inf.ethz.ch>
//
// coordinator   : ETH Zurich (Bernd Gaertner <gaertner@inf.ethz.ch>)
//
// Compute all furthest neighbors for the vertices of a convex polygon
// ============================================================================

#line 1408 "mon_search.aw"
#line 54 "code_formatting.awi"
#if ! (ALL_FURTHEST_NEIGHBORS_TRAITS_2_H)
#define ALL_FURTHEST_NEIGHBORS_TRAITS_2_H 1

#line 271 "afn.awi"
#ifndef CGAL_BASIC_H
#include <CGAL/basic.h>
#endif // CGAL_BASIC_H
#ifndef CGAL_OPTIMISATION_ASSERTIONS_H
#include <CGAL/optimisation_assertions.h>
#endif // CGAL_OPTIMISATION_ASSERTIONS_H
//!!! this should go into function_objects.h
#ifndef CGAL_PROTECT_FUNCTIONAL
#include <functional>
#define CGAL_PROTECT_FUNCTIONAL
#endif
#ifndef CGAL_SQUARED_DISTANCE_2_H
#include <CGAL/squared_distance_2.h>
#endif // CGAL_SQUARED_DISTANCE_2_H

#line 46 "code_formatting.awi"
CGAL_BEGIN_NAMESPACE
#line 283 "afn.awi"

template < class T1, class T2 >
struct Squared_distance
: public CGAL_STD::binary_function< T1, T2, typename T1::R::FT >
{
  typename T1::R::FT
  operator()( const T1& t1, const T2& t2) const
  { return squared_distance( t1, t2); }
};

#line 50 "code_formatting.awi"
CGAL_END_NAMESPACE
#line 294 "afn.awi"
#ifdef CGAL_CFG_NO_MEMBER_TEMPLATES
#ifndef CGAL_PROTECT_VECTOR
#include <vector>
#define CGAL_PROTECT_VECTOR
#endif
#endif
#line 46 "code_formatting.awi"
CGAL_BEGIN_NAMESPACE
#line 301 "afn.awi"

template < class _R >
class All_furthest_neighbors_traits {
public:
  typedef _R                                    R;
  typedef Point_2< R >                          Point_2;
  typedef Squared_distance< Point_2, Point_2 >  Distance;

#ifdef CGAL_CFG_NO_MEMBER_TEMPLATES
  typedef typename std::vector< Point_2 >::iterator
    RandomAccessIC;
  typedef typename std::vector< int >::reverse_iterator
    OutputIterator;
#endif

  #line 130 "traits.awi"
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
    typedef vector< Point_2 >            Cont;
    typedef Polygon_2< P_traits, Cont >  Polygon_2;
  
    Polygon_2 p( points_begin, points_end);
    return p.is_convex();
  } // is_convex( points_begin, points_end)
#line 317 "afn.awi"
};

template < class RandomAccessIC, class OutputIterator >
inline
OutputIterator
all_furthest_neighbors( RandomAccessIC points_begin,
                        RandomAccessIC points_end,
                        OutputIterator o)
{
  return
  _CGAL_all_furthest_neighbors(
    points_begin,
    points_end,
    o,
    std::value_type( points_begin));
} // all_furthest_neighbors( ... )

template < class RandomAccessIC,
           class OutputIterator,
           class R >
inline
OutputIterator
_CGAL_all_furthest_neighbors( RandomAccessIC points_begin,
                              RandomAccessIC points_end,
                              OutputIterator o,
                              Point_2< R >*)
{
  return
  all_furthest_neighbors(
    points_begin,
    points_end,
    o,
    All_furthest_neighbors_traits< R >());
} // _CGAL_all_furthest_neighbors( ... )

#line 50 "code_formatting.awi"
CGAL_END_NAMESPACE
#line 353 "afn.awi"

#endif // ! (ALL_FURTHEST_NEIGHBORS_TRAITS_2_H)

#line 12 "code_formatting.awi"
// ----------------------------------------------------------------------------
// ** EOF
// ----------------------------------------------------------------------------

