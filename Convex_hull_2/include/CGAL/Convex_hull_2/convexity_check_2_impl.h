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
// SPDX-License-Identifier: GPL-3.0+
// 
//
// Author(s)     : Stefan Schirra

#ifndef CGAL_CONVEXITY_CHECK_2_C
#define CGAL_CONVEXITY_CHECK_2_C

#include <CGAL/license/Convex_hull_2.h>


#include <CGAL/algorithm.h>
#include <algorithm>
#include <boost/bind.hpp>

namespace CGAL {

template <class ForwardIterator, class Traits>
bool
is_ccw_strongly_convex_2( ForwardIterator first, ForwardIterator last, 
                          const Traits& ch_traits)
{
  typedef  typename Traits::Less_xy_2      Less_xy;
  typedef  typename Traits::Left_turn_2    Left_turn;
  // added 
  typedef  typename Traits::Equal_2        Equal_2;
  
  Less_xy  smaller_xy    = ch_traits.less_xy_2_object();
  Left_turn left_turn    = ch_traits.left_turn_2_object();
  Equal_2  equal_points  = ch_traits.equal_2_object();   

  ForwardIterator iter1;
  ForwardIterator iter2;
  ForwardIterator iter3;

  if ( first == last) return true;

  iter2 = first;
  iter3 = ++iter2;

  if (iter3 == last ) return true;

  ++iter3;

  if (iter3 == last ) return (! equal_points(*first,*iter2) );

  iter1 = first;
  short int f = 0;

  while (iter3 != last) 
  {
      if ( !left_turn( *iter1, *iter2, *iter3 ) ) return false; 
      if ( smaller_xy( *iter2, *iter1 ) && smaller_xy( *iter2, *iter3 )) ++f;

      ++iter1;
      ++iter2;
      ++iter3;
  }

  iter3 = first;
  if ( !left_turn( *iter1, *iter2, *iter3 ) ) return false; 
  if ( smaller_xy( *iter2, *iter1 ) && smaller_xy( *iter2, *iter3 )) ++f;


  iter1 = iter2;
  iter2 = first;
  ++iter3;
  if ( !left_turn( *iter1, *iter2, *iter3 ) ) return false; 
  if ( smaller_xy( *iter2, *iter1 ) && smaller_xy( *iter2, *iter3 )) ++f;


  return ( f > 1 ) ? false : true;
}

template <class ForwardIterator, class Traits>
bool
is_cw_strongly_convex_2( ForwardIterator first, ForwardIterator last, 
                         const Traits& ch_traits)
{
  typedef  typename Traits::Less_xy_2      Less_xy;
  typedef  typename Traits::Equal_2        Equal_2;  

  Less_xy  smaller_xy    = ch_traits.less_xy_2_object();
  Equal_2  equal_points  = ch_traits.equal_2_object();  

  ForwardIterator iter1;
  ForwardIterator iter2;
  ForwardIterator iter3;

  if ( first == last) return true;

  iter2 = first;
  iter3 = ++iter2;

  if (iter3 == last ) return true;

  ++iter3;

  if (iter3 == last ) return (! equal_points(*first,*iter2) );

  iter1 = first;
  short int f = 0;

  while (iter3 != last) 
  {
      if ( !left_turn( *iter2, *iter1, *iter3 ) ) return false;
      if ( smaller_xy( *iter2, *iter1 ) && smaller_xy( *iter2, *iter3 )) ++f;

      ++iter1;
      ++iter2;
      ++iter3;
  }

  iter3 = first;
  if ( !left_turn( *iter2, *iter1, *iter3 ) ) return false;
  if ( smaller_xy( *iter2, *iter1 ) && smaller_xy( *iter2, *iter3 )) ++f;


  iter1 = iter2;
  iter2 = first;
  ++iter3;
  if ( !left_turn( *iter2, *iter1, *iter3 ) ) return false;
  if ( smaller_xy( *iter2, *iter1 ) && smaller_xy( *iter2, *iter3 )) ++f;


  return ( f > 1 ) ? false : true;
}

template <class ForwardIterator1, class ForwardIterator2, class Traits>
bool
ch_brute_force_check_2(ForwardIterator1 first1, ForwardIterator1 last1,
                            ForwardIterator2 first2, ForwardIterator2 last2,
                            const Traits&  ch_traits)
{
  using namespace boost;

  typedef    typename Traits::Left_turn_2    Left_of_line;
  ForwardIterator1 iter11;
  ForwardIterator2 iter21;
  ForwardIterator2 iter22;

  if ( first1 == last1) return true;

  if ( first2 == last2) return false;

  if ( cpp11::next(first2) == last2 )
  {
      while (first1 != last1)
      {
          if ( *first1++ != *first2 ) return false;
      }
      return true;
  }

  Left_of_line  left_turn = ch_traits.left_turn_2_object();
  iter22 = first2;
  iter21 = iter22++;
  while (iter22 != last2)
  {
      iter11 = std::find_if( first1, last1, 
                             bind(left_turn, *iter22++, *iter21++, _1) );
      if (iter11 != last1 ) return false;
  }

  iter11 = std::find_if( first1, last1, 
                         bind(left_turn, *first2, *iter21, _1) );
  if (iter11 != last1 ) return false;
  return true;
}

template <class ForwardIterator1, class ForwardIterator2, class Traits>
bool
ch_brute_force_chain_check_2(ForwardIterator1 first1, 
                                  ForwardIterator1 last1,
                                  ForwardIterator2 first2, 
                                  ForwardIterator2 last2,
                                  const Traits& ch_traits )
{
  using namespace boost;

  typedef  typename Traits::Left_turn_2     Left_turn_2;  
  
  ForwardIterator1 iter11;
  ForwardIterator2 iter21;
  ForwardIterator2 iter22;

  if ( first1 == last1) return true;

  if ( first2 == last2) return false;

  if ( cpp11::next(first2) == last2 ) return true;

  Left_turn_2  left_turn = ch_traits.left_turn_2_object();
  iter22 = first2;
  iter21 = iter22++;
  while (iter22 != last2)
  {
      iter11 = std::find_if( first1, last1, 
                             bind(left_turn, *iter22++, *iter21++, _1) );
      if (iter11 != last1 ) return false;
  }

  return true;
}

} //namespace CGAL

#endif // CGAL_CONVEXITY_CHECK_2_C
