// Copyright (c) 1999  Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
//
// Author(s)     : Stefan Schirra

#ifndef CGAL_CH_AKL_TOUSSAINT_C
#define CGAL_CH_AKL_TOUSSAINT_C

#ifndef CGAL_CH_NO_POSTCONDITIONS
#include <CGAL/convexity_check_2.h>
#endif // CGAL_CH_NO_POSTCONDITIONS

#include <CGAL/Convex_hull_2/ch_assertions.h>
#include <CGAL/ch_selected_extreme_points_2.h>
#include <CGAL/ch_graham_andrew.h>
#include <CGAL/algorithm.h>
#include <CGAL/IO/Tee_for_output_iterator.h>
#include <boost/bind.hpp>
#include <CGAL/tuple.h>
#include <CGAL/utility.h>
#include <iterator>

namespace CGAL {

  
namespace internal{
  
template <class ForwardIterator, class Traits>
inline
cpp11::tuple<ForwardIterator,ForwardIterator,ForwardIterator,ForwardIterator>
ch_nswe_point_with_order( ForwardIterator first, ForwardIterator last,
                          ForwardIterator& n,
                          ForwardIterator& s,
                          ForwardIterator& w,
                          ForwardIterator& e,
                          const Traits& ch_traits,
                          std::forward_iterator_tag  )
{
  typename Traits::Less_xy_2    
      lexicographically_xy_smaller = ch_traits.less_xy_2_object();
  typename Traits::Less_yx_2    
      lexicographically_yx_smaller = ch_traits.less_yx_2_object();
  n = s = w = e = first;
  unsigned i=0;
  //array use to track the position of w,e,n,s in the range. first is for the position, second to track
  //the original position after sorting
  std::pair<unsigned,unsigned> positions[4]={
    std::make_pair(0,0),
    std::make_pair(0,1),
    std::make_pair(0,2),
    std::make_pair(0,3) };

  while ( first != last )
  {
      if ( lexicographically_xy_smaller( *first, *w ))  { w = first; positions[0].first=i; }
      if ( lexicographically_xy_smaller( *e, *first ))  { e = first; positions[1].first=i; }
      if ( lexicographically_yx_smaller( *n, *first ))  { n = first; positions[2].first=i; }
      if ( lexicographically_yx_smaller( *first, *s ))  { s = first; positions[3].first=i; }
      ++first;
      ++i;
  }
  ForwardIterator iterators[4]={w,e,n,s};
  std::sort(positions,positions+4);
  
  return cpp11::make_tuple( 
    iterators[positions[0].second],
    iterators[positions[1].second],
    iterators[positions[2].second],
    iterators[positions[3].second]
  );
}
  
  
template <class RandomAccessIterator, class Traits>
inline
cpp11::tuple<RandomAccessIterator,RandomAccessIterator,RandomAccessIterator,RandomAccessIterator>
ch_nswe_point_with_order( RandomAccessIterator first, RandomAccessIterator last,
                          RandomAccessIterator& n,
                          RandomAccessIterator& s,
                          RandomAccessIterator& w,
                          RandomAccessIterator& e,
                          const Traits& ch_traits,
                          std::random_access_iterator_tag)
{
  ch_nswe_point(first,last,n,s,w,e,ch_traits);
  RandomAccessIterator iterators[4]={w,e,n,s};
  std::sort(iterators,iterators+4);
  
  return cpp11::make_tuple( 
    iterators[0],
    iterators[1],
    iterators[2],
    iterators[3]
  );
}

//this function does the same as ch_nswe_point but return the iterators n,s,w,e
//sorted according to their positions in the range [first,last]
template <class ForwardIterator, class Traits>
inline
cpp11::tuple<ForwardIterator,ForwardIterator,ForwardIterator,ForwardIterator>
ch_nswe_point_with_order( ForwardIterator first, ForwardIterator last,
                          ForwardIterator& n,
                          ForwardIterator& s,
                          ForwardIterator& w,
                          ForwardIterator& e,
                          const Traits& ch_traits)
{
  return ch_nswe_point_with_order(first,last,n,s,w,e,ch_traits,typename std::iterator_traits<ForwardIterator>::iterator_category());
}

template <class ForwardIterator,class Traits>
inline
void ch_akl_toussaint_assign_points_to_regions(ForwardIterator first, ForwardIterator last, 
                                               const typename Traits::Left_turn_2&  left_turn,
                                               ForwardIterator e,
                                               ForwardIterator w,
                                               ForwardIterator n,
                                               ForwardIterator s,
                                               std::vector< typename Traits::Point_2 >& region1,
                                               std::vector< typename Traits::Point_2 >& region2,
                                               std::vector< typename Traits::Point_2 >& region3,
                                               std::vector< typename Traits::Point_2 >& region4,
                                               const Traits&)
{
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
}

template <class ForwardIterator,class Traits>
inline
void ch_akl_toussaint_assign_points_to_regions_deg(ForwardIterator first, ForwardIterator last, 
                                                   const typename Traits::Left_turn_2&  left_turn,
                                                   ForwardIterator e,
                                                   ForwardIterator w,
                                                   ForwardIterator n,
                                                   ForwardIterator s,
                                                   std::vector< typename Traits::Point_2 >& region1,
                                                   std::vector< typename Traits::Point_2 >& region2,
                                                   std::vector< typename Traits::Point_2 >& region3,
                                                   std::vector< typename Traits::Point_2 >& region4,
                                                   int duplicated_exteme_points,
                                                   const Traits& traits)
{
  std::vector< typename Traits::Point_2 >& r1 = (s==w?region2:region1);
  std::vector< typename Traits::Point_2 >& r3 = (n==e?region4:region3);
  switch(duplicated_exteme_points){
    case 2:
    {
      typename Traits::Orientation_2 orient = traits.orientation_2_object();
      for ( ; first != last; ++first )
      {   
        switch( orient(*e,*w,*first) ){
          case LEFT_TURN: 
            r1.push_back( *first );
          break;
          case RIGHT_TURN: 
            r3.push_back( *first );
          break;
          default:
          break;
        }
      }
      break;
    }
    default: //this is case 1
    if (s==w || s==e){
      for ( ; first != last; ++first )
      {
        if ( left_turn(*e, *w, *first ) )
          r1.push_back( *first );
        else
        {
            if ( left_turn( *n, *e, *first ) )       region3.push_back( *first );
            else if ( left_turn( *w, *n, *first ) )  region4.push_back( *first );
        }
      }
    }
    else{
      for ( ; first != last; ++first )
      {
        //note that e!=w and s!=n except if the convex hull is a point (they are lexicographically sorted)
        if ( left_turn(*e, *w, *first ) )   
        {
            if (s!=w && left_turn( *s, *w, *first ) )       region1.push_back( *first );
            else if (e!=s && left_turn( *e, *s, *first ) )  region2.push_back( *first );
        }
        else
          r3.push_back( *first );
      }
    }
  }
}

}//namespace internal

template <class ForwardIterator, class OutputIterator, class Traits>
OutputIterator
ch_akl_toussaint(ForwardIterator first, ForwardIterator last, 
                 OutputIterator  result,
                 const Traits&   ch_traits)
{
  using namespace boost;

  typedef  typename Traits::Point_2                    Point_2;    
  typedef  typename Traits::Left_turn_2                Left_of_line;
  // added 
  typedef  typename Traits::Equal_2                    Equal_2;
  
  Left_of_line  left_turn        = ch_traits.left_turn_2_object();
  Equal_2       equal_points     = ch_traits.equal_2_object();  

  if (first == last) return result;
  ForwardIterator n, s, e, w;
  cpp11::tuple<ForwardIterator,ForwardIterator,ForwardIterator,ForwardIterator> ranges=
    internal::ch_nswe_point_with_order( first, last, n, s, w, e, ch_traits);
  
  if (equal_points(*n, *s) )
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

  CGAL_ch_postcondition_code( ForwardIterator save_first = first; )
  
  int duplicated_exteme_points =  (cpp11::get<0>(ranges)==cpp11::get<1>(ranges)?1:0) +
                                  (cpp11::get<1>(ranges)==cpp11::get<2>(ranges)?1:0) +
                                  (cpp11::get<2>(ranges)==cpp11::get<3>(ranges)?1:0);
  
  //several calls to avoid filter failures when using n,s,e,w
  if (duplicated_exteme_points)
  {
    internal::ch_akl_toussaint_assign_points_to_regions_deg(first,cpp11::get<0>(ranges),left_turn,e,w,n,s,region1,region2,region3,region4,duplicated_exteme_points,ch_traits);
    
    if ( cpp11::get<0>(ranges)!=cpp11::get<1>(ranges) )
      internal::ch_akl_toussaint_assign_points_to_regions_deg(cpp11::next(cpp11::get<0>(ranges)),cpp11::get<1>(ranges),left_turn,e,w,n,s,region1,region2,region3,region4,duplicated_exteme_points,ch_traits);
    
    if ( cpp11::get<1>(ranges)!=cpp11::get<2>(ranges) )
      internal::ch_akl_toussaint_assign_points_to_regions_deg(cpp11::next(cpp11::get<1>(ranges)),cpp11::get<2>(ranges),left_turn,e,w,n,s,region1,region2,region3,region4,duplicated_exteme_points,ch_traits);
    
    if ( cpp11::get<2>(ranges)!=cpp11::get<3>(ranges) )
      internal::ch_akl_toussaint_assign_points_to_regions_deg(cpp11::next(cpp11::get<2>(ranges)),cpp11::get<3>(ranges),left_turn,e,w,n,s,region1,region2,region3,region4,duplicated_exteme_points,ch_traits);
    
    internal::ch_akl_toussaint_assign_points_to_regions_deg(cpp11::next(cpp11::get<3>(ranges)),last,left_turn,e,w,n,s,region1,region2,region3,region4,duplicated_exteme_points,ch_traits);
  }
  else{
    internal::ch_akl_toussaint_assign_points_to_regions(first,cpp11::get<0>(ranges),left_turn,e,w,n,s,region1,region2,region3,region4,ch_traits);
    internal::ch_akl_toussaint_assign_points_to_regions(cpp11::next(cpp11::get<0>(ranges)),cpp11::get<1>(ranges),left_turn,e,w,n,s,region1,region2,region3,region4,ch_traits);
    internal::ch_akl_toussaint_assign_points_to_regions(cpp11::next(cpp11::get<1>(ranges)),cpp11::get<2>(ranges),left_turn,e,w,n,s,region1,region2,region3,region4,ch_traits);
    internal::ch_akl_toussaint_assign_points_to_regions(cpp11::next(cpp11::get<2>(ranges)),cpp11::get<3>(ranges),left_turn,e,w,n,s,region1,region2,region3,region4,ch_traits);
    internal::ch_akl_toussaint_assign_points_to_regions(cpp11::next(cpp11::get<3>(ranges)),last,left_turn,e,w,n,s,region1,region2,region3,region4,ch_traits);
  }
  
  #if defined(CGAL_CH_NO_POSTCONDITIONS) || defined(CGAL_NO_POSTCONDITIONS) \
    || defined(NDEBUG)
  OutputIterator  res(result);
  #else
  Tee_for_output_iterator<OutputIterator,Point_2> res(result);
  #endif // no postconditions ...
  std::sort( cpp11::next(region1.begin() ), region1.end(), 
             ch_traits.less_xy_2_object() );
  std::sort( cpp11::next(region2.begin() ), region2.end(), 
             ch_traits.less_xy_2_object() );
  std::sort( cpp11::next(region3.begin() ), region3.end(),
             boost::bind(ch_traits.less_xy_2_object(), _2, _1) );
  std::sort( cpp11::next(region4.begin() ), region4.end(), 
             boost::bind(ch_traits.less_xy_2_object(), _2, _1) );

  if (! equal_points(*w,*s) )
  {
      region1.push_back( *s );
      ch__ref_graham_andrew_scan( region1.begin(), region1.end(), 
                                       res, ch_traits);
  }
  if (! equal_points(*s,*e) )
  {
      region2.push_back( *e );
      ch__ref_graham_andrew_scan( region2.begin(), region2.end(),
                                       res, ch_traits);
  }
  if (! equal_points(*e,*n) )
  {
      region3.push_back( *n );
      ch__ref_graham_andrew_scan( region3.begin(), region3.end(),
                                       res, ch_traits);
  }
  if (! equal_points(*n,*w) )
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

} //namespace CGAL

#endif // CGAL_CH_AKL_TOUSSAINT_C
