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
// file          : ch_bykat.C
// package       : Convex_hull (3.3)
// maintainer    : Stefan Schirra <stschirr@mpi-sb.mpg.de>
// source        : convex_hull_2.lw
// revision      : 3.3
// revision_date : 03 Aug 2000
// author(s)     : Stefan Schirra
//
// coordinator   : MPI, Saarbruecken
// ======================================================================


#ifndef CGAL_CH_BYKAT_C
#define CGAL_CH_BYKAT_C

#ifndef CGAL_CH_BYKAT_H
#include <CGAL/ch_eddy.h>
#endif // CGAL_CH_BYKAT_H
#include <CGAL/functional.h>

CGAL_BEGIN_NAMESPACE
template <class InputIterator, class OutputIterator, class Traits>
OutputIterator
ch_bykat(InputIterator first, InputIterator last, 
              OutputIterator  result,
              const Traits& ch_traits)
{
  typedef typename Traits::Point_2                         Point_2;
  typedef typename Traits::Leftturn_2                      Leftturn_2;
  typedef typename Traits::Less_signed_distance_to_line_2  Less_dist;

  if (first == last) return result;

  std::vector< Point_2 >       P;      // Points in subsets
  std::vector< Point_2 >       H;      // right endpoints of subproblems
  P.reserve(16);
  H.reserve(16);
  
  typedef typename std::vector< Point_2 >::iterator   PointIterator;
  std::vector< PointIterator > L;      // start of subset range
  std::vector< PointIterator > R;      // end of subset range
  L.reserve(16);
  R.reserve(16);
  PointIterator           l;
  PointIterator           r;
  Point_2                 a,b,c;
  
  std::copy(first,last,std::back_inserter(P));
  ch_we_point(P.begin(), P.end(), l, r, ch_traits);
  a = *l;
  b = *r;
  if ( a == b) 
  {
      *result = a;  ++result;
      return result;
  }
  #if defined(CGAL_CH_NO_POSTCONDITIONS) || defined(CGAL_NO_POSTCONDITIONS) \
    || defined(NDEBUG)
  OutputIterator  res(result);
  #else
  Tee_for_output_iterator<OutputIterator,Point_2> res(result);
  #endif // no postconditions ...
  H.push_back( a );
  L.push_back( P.begin() );
  Leftturn_2 left_turn = ch_traits.leftturn_2_object();
  R.push_back( l = std::partition( P.begin(), P.end(), 
                                   bind_1(bind_1(left_turn,a),b) ) );
  r = std::partition( l, P.end(), bind_1(bind_1(left_turn,b),a) );
  
  Less_dist less_dist = ch_traits.less_signed_distance_to_line_2_object();
  for (;;)
  {
      if ( l != r)
      {
          c = *std::min_element( l, r, bind_1(bind_1(less_dist, a), b));
          H.push_back( b );
          L.push_back( l );
          R.push_back( l = std::partition(l, r, 
                                          bind_1(bind_1(left_turn,b),c)));
          r = std::partition(l, r, bind_1(bind_1(left_turn,c),a));
          b = c; 
      }
      else
      {
          *res = a;  ++res;
          if ( L.empty() ) break;
          a = b;
          b = H.back(); H.pop_back();
          l = L.back(); L.pop_back();
          r = R.back(); R.pop_back();
      }
  }
  CGAL_ch_postcondition( \
      is_ccw_strongly_convex_2( res.output_so_far_begin(), \
                                     res.output_so_far_end(), \
                                     ch_traits));
  CGAL_ch_expensive_postcondition( \
      ch_brute_force_check_2( \
          P.begin(), P.end(), \
          res.output_so_far_begin(), res.output_so_far_end(), \
          ch_traits));
  #if defined(CGAL_CH_NO_POSTCONDITIONS) || defined(CGAL_NO_POSTCONDITIONS) \
    || defined(NDEBUG)
  return res;
  #else
  return res.to_output_iterator();
  #endif // no postconditions ...
}

#define CGAL_ch_THRESHOLD 10
template <class InputIterator, class OutputIterator, class Traits>
OutputIterator
ch_bykat_with_threshold(InputIterator   first, InputIterator last, 
                             OutputIterator  result,
                             const Traits&   ch_traits)
{
  typedef typename Traits::Point_2               Point_2;
  typedef typename Traits::Leftturn_2            Leftturn_2;
  typedef typename Traits::Less_signed_distance_to_line_2     
                                                 Less_dist;
  typedef typename Traits::Less_xy_2             Less_xy;
  typedef ch_Binary_predicate_reversor< Point_2, Less_xy >
                                                 Greater_xy;
  typedef typename std::vector< Point_2 >::iterator   
                                                 PointIterator;

  if (first == last) return result;

  std::vector< Point_2 >       P;      // points in subsets
  std::vector< Point_2 >       H;      // right endpoints of subproblems
  P.reserve(16);
  H.reserve(16);
  std::vector< PointIterator > L;      // start of subset range
  std::vector< PointIterator > R;      // end of subset range
  L.reserve(16);
  R.reserve(16);
  PointIterator           l;
  PointIterator           r;
  Point_2                 a,b,c;
  PointIterator           Pbegin, Pend;
  
  P.push_back(Point_2() );
  std::copy(first,last,std::back_inserter(P));
  P.push_back(Point_2() );
  Pbegin = successor(P.begin());
  Pend   = predecessor(P.end());
  ch_we_point(Pbegin, Pend, l, r, ch_traits);
  a = *l;
  b = *r;
  if ( a == b) 
  {
      *result = a;  ++result;
      return result;
  }
  #if defined(CGAL_CH_NO_POSTCONDITIONS) || defined(CGAL_NO_POSTCONDITIONS) \
    || defined(NDEBUG)
  OutputIterator  res(result);
  #else
  Tee_for_output_iterator<OutputIterator,Point_2> res(result);
  #endif // no postconditions ...
  H.push_back( a );
  L.push_back( Pbegin );
  Leftturn_2 left_turn = ch_traits.leftturn_2_object();
  R.push_back( l = std::partition( Pbegin, Pend, 
                                   bind_1(bind_1(left_turn, a), b) ) );
  r = std::partition( l, Pend, bind_1(bind_1(left_turn,b),a) );
  
  Less_dist less_dist = ch_traits.less_signed_distance_to_line_2_object();
  for (;;)
  {
      if ( l != r)
      {
          if ( r-l > CGAL_ch_THRESHOLD )
          {
              c = *std::min_element( l, r, bind_1(bind_1(less_dist, a), b));
              H.push_back( b );
              L.push_back( l );
              R.push_back( l = std::partition(l, r, 
                           bind_1(bind_1(left_turn, b), c)) );
              r = std::partition(l, r, bind_1(bind_1(left_turn, c), a));
              b = c; 
          }
          else
          {
              std::swap( a, *--l);
              std::swap( b, *++r);
              if ( ch_traits.less_xy_2_object()(*l,*r) )
              {
                  std::sort(successor(l), r, 
                            ch_traits.less_xy_2_object() );
              }
              else
              {
                  std::sort(successor(l), r, 
                            Greater_xy(ch_traits.less_xy_2_object()) );
              }
              ch__ref_graham_andrew_scan(l, successor(r), res, ch_traits);
              std::swap( a, *l);
              std::swap( b, *r);
              if ( L.empty() ) break;
              a = b;
              b = H.back(); H.pop_back();
              l = L.back(); L.pop_back();
              r = R.back(); R.pop_back();
          }
              
      }
      else
      {
          *res = a;  ++res;
          if ( L.empty() ) break;
          a = b;
          b = H.back(); H.pop_back();
          l = L.back(); L.pop_back();
          r = R.back(); R.pop_back();
      }
  }
  CGAL_ch_postcondition( \
      is_ccw_strongly_convex_2( res.output_so_far_begin(), \
                                     res.output_so_far_end(), \
                                     ch_traits));
  CGAL_ch_expensive_postcondition( \
      ch_brute_force_check_2( \
          Pbegin, Pend, \
          res.output_so_far_begin(), res.output_so_far_end(), \
          ch_traits));
  #if defined(CGAL_CH_NO_POSTCONDITIONS) || defined(CGAL_NO_POSTCONDITIONS) \
    || defined(NDEBUG)
  return res;
  #else
  return res.to_output_iterator();
  #endif // no postconditions ...
}

CGAL_END_NAMESPACE

#endif // CGAL_CH_BYKAT_C


