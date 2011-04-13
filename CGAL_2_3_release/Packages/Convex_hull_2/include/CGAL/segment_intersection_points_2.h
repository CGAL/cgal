// ======================================================================
//
// Copyright (c) 1998 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
// release       : 
// release_date  : 
//
// file          : include/CGAL/segment_intersection_points_2.h
// package       : Convex_hull_2
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Stefan Schirra 
//
// coordinator   : MPI, Saarbruecken
// ======================================================================
#ifndef SEGMENT_INTERSECTION_POINTS_2_H
#define SEGMENT_INTERSECTION_POINTS_2_H

#include <CGAL/stl_extensions.h>

namespace CGAL {
/*
#include <CGAL/Segment_2_Segment_2_intersection.h>
template <class ForwardIterator, class OutputIterator, class R>
OutputIterator
si_brute_force(ForwardIterator first, ForwardIterator last,
                    OutputIterator result,
                    const R& )
{
  typedef Point_2<R>         Point;
  ForwardIterator inner, outer;
  Object i_obj;
  Point p;
  for ( outer = first; outer != last; ++outer)
      for ( inner = successor(outer); inner != last; ++inner)
      {
          i_obj = intersection( *outer, *inner);
          if ( assign( p, i_obj) )
              result++ = p;
      }
  return result;
}
*/

template <class R>
Point_2<R>
gp_linear_intersection(const Line_2<R>& l1,
                            const Line_2<R>& l2)
{
  return Point_2<R>( l1.b()*l2.c() - l2.b()*l1.c(),
                     l2.a()*l1.c() - l1.a()*l2.c(),
                     l1.a()*l2.b() - l2.a()*l1.b() );
}

template <class ForwardIterator, class OutputIterator, class R>
OutputIterator
si_brute_force_II(ForwardIterator first, ForwardIterator last,
                  OutputIterator result,
                  const R& )
{
  typedef Point_2<R>         Point;
  ForwardIterator inner, outer;
  Point p;
  for ( outer = first; outer != last; ++outer)
      for ( inner = successor(outer); inner != last; ++inner)
      {
          Point s1 = (*outer).source();
          Point e1 = (*outer).target();
          Point s2 = (*inner).source();
          Point e2 = (*inner).target();
          if ( (orientation( s1, e1, s2) != orientation( s1, e1, e2))
             &&(orientation( s2, e2, s1) != orientation( s2, e2, e1)))
          {
              result++ = gp_linear_intersection( Line_2<R>(s1,e1),
                                                      Line_2<R>(s2,e2));
          }
      }
  return result;
}

template <class ForwardIterator, class OutputIterator, class R>
OutputIterator
segment_intersection_points_2(ForwardIterator first, ForwardIterator last,
                                   OutputIterator result,
                                   const R& r)
{ return si_brute_force_II( first, last, result, r); }

} // namespace CGAL

#endif // SEGMENT_INTERSECTION_POINTS_2_H
