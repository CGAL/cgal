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
// file          : ch_akl_toussaint.C
// package       : Convex_hull (3.3)
// maintainer    : Stefan Schirra <stschirr@mpi-sb.mpg.de>
// source        : convex_hull_2.lw
// revision      : 3.3
// revision_date : 03 Aug 2000
// author(s)     : Stefan Schirra
//
// coordinator   : MPI, Saarbruecken
// ======================================================================


#ifndef CGAL_CH_AKL_TOUSSAINT_C
#define CGAL_CH_AKL_TOUSSAINT_C

#ifndef CGAL_CH_AKL_TOUSSAINT_H
#include <CGAL/ch_akl_toussaint.h>
#endif // CGAL_CH_AKL_TOUSSAINT_H

CGAL_BEGIN_NAMESPACE
template <class ForwardIterator, class OutputIterator, class Traits>
OutputIterator
ch_akl_toussaint(ForwardIterator first, ForwardIterator last, 
                      OutputIterator  result,
                      const Traits&   ch_traits)
{
  typedef  typename Traits::Point_2                    Point_2;    
  typedef  typename Traits::Leftturn_2                 Left_of_line;
  typedef  typename Traits::Less_xy_2                  Less_xy;
  typedef  ch_Binary_predicate_reversor< Point_2, Less_xy>
                                                       Greater_xy;
  typedef  typename Traits::Less_yx_2                  Less_yx;
  typedef  ch_Binary_predicate_reversor< Point_2, Less_yx>
                                                       Greater_yx;

  if (first == last) return result;
  ForwardIterator n, s, e, w;
  ch_nswe_point( first, last, n, s, w, e, ch_traits);
  if ( *n == *s )
  {
      *result = *w;  ++result;
      return result;
  }


  std::vector< Point_2 > region1;
  std::vector< Point_2 > region2;
  std::vector< Point_2 > region3;
  std::vector< Point_2 > region4;
  region1.reserve(16);
  region2.reserve(16);
  region3.reserve(16);
  region4.reserve(16);
  region1.push_back( *w);
  region2.push_back( *s);
  region3.push_back( *e);
  region4.push_back( *n);

  Left_of_line  left_turn = ch_traits.leftturn_2_object();

  CGAL_ch_postcondition_code( ForwardIterator save_first = first; )

  for ( ; first != last; ++first )
  {
      if ( left_turn(*e, *w, *first ) )   
      {
          if ( left_turn( *s, *w, *first ) )       region1.push_back( *first );
          else if ( left_turn( *e, *s, *first ) )  region2.push_back( *first );
      }
      else
      {
          if ( left_turn( *n, *e, *first ) )       region3.push_back( *first );
          else if ( left_turn( *w, *n, *first ) )  region4.push_back( *first );
      }
  }

  #if defined(CGAL_CH_NO_POSTCONDITIONS) || defined(CGAL_NO_POSTCONDITIONS) \
    || defined(NDEBUG)
  OutputIterator  res(result);
  #else
  Tee_for_output_iterator<OutputIterator,Point_2> res(result);
  #endif // no postconditions ...
  std::sort( successor(region1.begin() ), region1.end(), 
             ch_traits.less_xy_2_object() );
  std::sort( successor(region2.begin() ), region2.end(), 
             ch_traits.less_xy_2_object() );
  std::sort( successor(region3.begin() ), region3.end(), 
             Greater_xy(ch_traits.less_xy_2_object()) );
  std::sort( successor(region4.begin() ), region4.end(), 
             Greater_xy(ch_traits.less_xy_2_object()) );

  if ( *w != *s )
  {
      region1.push_back( *s );
      ch__ref_graham_andrew_scan( region1.begin(), region1.end(), 
                                       res, ch_traits);
  }
  if ( *s != *e )
  {
      region2.push_back( *e );
      ch__ref_graham_andrew_scan( region2.begin(), region2.end(),
                                       res, ch_traits);
  }
  if ( *e != *n )
  {
      region3.push_back( *n );
      ch__ref_graham_andrew_scan( region3.begin(), region3.end(),
                                       res, ch_traits);
  }
  if ( *n != *w )
  {
      region4.push_back( *w );
      ch__ref_graham_andrew_scan( region4.begin(), region4.end(),
                                       res, ch_traits);
  }

  CGAL_ch_postcondition_code( first = save_first; )
  CGAL_ch_postcondition( \
      is_ccw_strongly_convex_2( res.output_so_far_begin(), \
                                     res.output_so_far_end(), \
                                     ch_traits));
  CGAL_ch_expensive_postcondition( \
      ch_brute_force_check_2( \
          first, last, \
          res.output_so_far_begin(), res.output_so_far_end(), \
          ch_traits)
  );
  #if defined(CGAL_CH_NO_POSTCONDITIONS) || defined(CGAL_NO_POSTCONDITIONS) \
    || defined(NDEBUG)
  return res;
  #else
  return res.to_output_iterator();
  #endif // no postconditions ...

}

CGAL_END_NAMESPACE

#endif // CGAL_CH_AKL_TOUSSAINT_C


