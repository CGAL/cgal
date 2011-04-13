// ======================================================================
//
// Copyright (c) 1999 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
// release       : 
// release_date  : 2000, August 03
//
// file          : ch_eddy.C
// package       : Convex_hull (3.3)
// maintainer    : Stefan Schirra <stschirr@mpi-sb.mpg.de>
// source        : convex_hull_2.lw
// revision      : 3.3
// revision_date : 03 Aug 2000
// author(s)     : Stefan Schirra
//
// coordinator   : MPI, Saarbruecken
// ======================================================================


#ifndef CGAL_CH_EDDY_C
#define CGAL_CH_EDDY_C

#ifndef CGAL_CH_EDDY_H
#include <CGAL/ch_eddy.h>
#endif // CGAL_CH_EDDY_H
#include <CGAL/functional.h>

CGAL_BEGIN_NAMESPACE
template <class List, class ListIterator, class Traits>
void
ch__recursive_eddy(List& L, 
                        ListIterator  a_it, ListIterator  b_it, 
                        const Traits& ch_traits)
{
  typedef  typename Traits::Point_2                         Point_2;    
  typedef  typename Traits::Leftturn_2                      Leftturn_2;
  typedef  typename Traits::Less_signed_distance_to_line_2  Less_dist;

  Leftturn_2 left_turn = ch_traits.leftturn_2_object();
  CGAL_ch_precondition( \
    std::find_if(a_it, b_it, \
            bind_1(bind_1(left_turn, *b_it), *a_it)) \
    != b_it );


  ListIterator f_it = successor(a_it);
  Less_dist less_dist = ch_traits.less_signed_distance_to_line_2_object();
  ListIterator 
      c_it = std::min_element( f_it, b_it,  // max before
                               bind_1(bind_1(less_dist, *a_it), *b_it));
  Point_2 c = *c_it;

  c_it = std::partition(f_it, b_it, bind_1(bind_1(left_turn, c), *a_it));
  f_it = std::partition(c_it, b_it, bind_1(bind_1(left_turn, *b_it), c));
  c_it = L.insert(c_it, c);
  L.erase( f_it, b_it );

  if ( successor(a_it) != c_it )
  {
      ch__recursive_eddy( L, a_it, c_it, ch_traits);
  }
  if ( successor(c_it) != b_it )
  {
      ch__recursive_eddy( L, c_it, b_it, ch_traits);
  }


}

template <class InputIterator, class OutputIterator, class Traits>
OutputIterator
ch_eddy(InputIterator first, InputIterator last,
             OutputIterator  result,
             const Traits& ch_traits)
{
  typedef  typename Traits::Point_2                         Point_2;    
  typedef  typename Traits::Leftturn_2                      Leftturn_2;

  if (first == last) return result;
  std::list< Point_2 >   L;
  std::copy( first, last, std::back_inserter(L) );

  typedef typename std::list< Point_2 >::iterator  list_iterator;
  list_iterator   w, e;
  ch_we_point(L.begin(), L.end(), w, e, ch_traits);
  Point_2 wp = *w;
  Point_2 ep = *e;
  if ( wp == ep )
  {
      *result = wp;  ++result;
      return result;
  }

  L.erase(w);
  L.erase(e);
  Leftturn_2 left_turn = ch_traits.leftturn_2_object();
  e = std::partition(L.begin(), L.end(), 
                     bind_1(bind_1(left_turn, ep), wp) );
  L.push_front(wp);
  e = L.insert(e, ep);

  if ( successor(L.begin()) != e )
  {
      ch__recursive_eddy( L, L.begin(), e, ch_traits);
  }
  w = std::find_if( e, L.end(), bind_1(bind_1(left_turn, wp), ep) );
  if ( w == L.end() )
  {
      L.erase( ++e, L.end() );
      return std::copy( L.begin(), L.end(), result );
  }
  w = L.insert(L.end(), wp);
  ch__recursive_eddy( L, e, w, ch_traits);


  CGAL_ch_postcondition( \
      is_ccw_strongly_convex_2( L.begin(), w, ch_traits) );
  CGAL_ch_expensive_postcondition( \
      ch_brute_force_check_2( first, last, \
                                   L.begin(), w, ch_traits ) );

  return std::copy( L.begin(), w, result );

}

CGAL_END_NAMESPACE

#endif // CGAL_CH_EDDY_C

