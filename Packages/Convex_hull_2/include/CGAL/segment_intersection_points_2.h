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

#include <CGAL/basic.h>
#include <iterator>

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

template <class ForwardIterator, class OutputIterator, class Traits>
OutputIterator
si_brute_force_II(ForwardIterator first, ForwardIterator last,
                  OutputIterator result,
                  const Traits& traits)
{
  typedef typename Traits::Point_2         Point;
  typedef typename Traits::Line_2          Line;
  typedef typename Traits::Orientation_2   Orientation;
  Orientation orientation = traits.orientation_2_object();

  for ( ForwardIterator outer = first; outer != last; ++outer)
      for ( ForwardIterator inner = successor(outer); inner != last; ++inner)
      {
          Point s1 = (*outer).source();
          Point e1 = (*outer).target();
          Point s2 = (*inner).source();
          Point e2 = (*inner).target();
          if ( (orientation( s1, e1, s2) != orientation( s1, e1, e2))
             &&(orientation( s2, e2, s1) != orientation( s2, e2, e1)))
          {
              Line l1( s1, e1);
              Line l2( s2, e2);
              result++ =  Point( l1.b()*l2.c() - l2.b()*l1.c(),
                                 l2.a()*l1.c() - l1.a()*l2.c(),
                                 l1.a()*l2.b() - l2.a()*l1.b());
          }
      }
  return result;
}

template <class ForwardIterator, class OutputIterator, class Traits>
OutputIterator
segment_intersection_points_2(ForwardIterator first, ForwardIterator last,
                              OutputIterator result,
                              const Traits& traits)
{ return si_brute_force_II( first, last, result, traits); }

template <class ForwardIterator, class OutputIterator, class Traits>
OutputIterator
segment_intersection_points_2(ForwardIterator first, ForwardIterator last,
                              OutputIterator result)
{
    typedef std::iterator_traits<ForwardIterator> ITraits;
    typedef typename ITraits::value_type          value_type;
    typedef CGAL::Kernel_traits<value_type>       KTraits;
    typedef typename KTraits::Kernel              Kernel;
    return segment_intersection_points_2( first, last, result, Kernel()); 
}

} // namespace CGAL

#endif // SEGMENT_INTERSECTION_POINTS_2_H
