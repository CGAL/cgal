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

#if ! (CGAL_ALL_FURTHEST_NEIGHBORS_TRAITS_2_H)
#define CGAL_ALL_FURTHEST_NEIGHBORS_TRAITS_2_H 1

#ifndef CGAL_OPTIMISATION_ASSERTIONS_H
#include <CGAL/optimisation_assertions.h>
#endif // CGAL_OPTIMISATION_ASSERTIONS_H
//!!! this should go into function_objects.h
#ifndef CGAL_PROTECT_FUNCTION_H
#include <function.h>
#define CGAL_PROTECT_FUNCTION_H
#endif // CGAL_PROTECT_FUNCTION_H
#ifndef CGAL_SQUARED_DISTANCE_2_H
#include <CGAL/squared_distance_2.h>
#endif // CGAL_SQUARED_DISTANCE_2_H
template < class T1, class T2 >
struct CGAL_Squared_distance
: public binary_function< T1, T2, typename T1::R::FT >
{
  typename T1::R::FT
  operator()( const T1& t1, const T2& t2) const
  { return CGAL_squared_distance( t1, t2); }
};

#ifdef CGAL_CFG_NO_MEMBER_TEMPLATES
#ifndef CGAL_PROTECT_VECTOR_H
#include <vector.h>
#define CGAL_PROTECT_VECTOR_H
#endif // CGAL_PROTECT_VECTOR_H
#endif

template < class _R >
class CGAL_All_furthest_neighbors_traits {
public:
  typedef _R                                         R;
  typedef CGAL_Point_2< R >                          Point_2;
  typedef CGAL_Squared_distance< Point_2, Point_2 >  Distance;

#ifdef CGAL_CFG_NO_MEMBER_TEMPLATES
  typedef typename vector< Point_2 >::iterator
    RandomAccessIC;
  typedef typename vector< int >::reverse_iterator
    OutputIterator;
#endif

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

template < class RandomAccessIC, class OutputIterator >
inline
OutputIterator
CGAL_all_furthest_neighbors( RandomAccessIC points_begin,
                             RandomAccessIC points_end,
                             OutputIterator o)
{
  return
  _CGAL_all_furthest_neighbors(
    points_begin,
    points_end,
    o,
    value_type( points_begin));
} // CGAL_all_furthest_neighbors( ... )

template < class RandomAccessIC,
           class OutputIterator,
           class R >
inline
OutputIterator
_CGAL_all_furthest_neighbors( RandomAccessIC points_begin,
                              RandomAccessIC points_end,
                              OutputIterator o,
                              CGAL_Point_2< R >*)
{
  return
  CGAL_all_furthest_neighbors(
    points_begin,
    points_end,
    o,
    CGAL_All_furthest_neighbors_traits< R >());
} // _CGAL_all_furthest_neighbors( ... )

#endif // ! (CGAL_ALL_FURTHEST_NEIGHBORS_TRAITS_2_H)

// ----------------------------------------------------------------------------
// ** EOF
// ----------------------------------------------------------------------------

