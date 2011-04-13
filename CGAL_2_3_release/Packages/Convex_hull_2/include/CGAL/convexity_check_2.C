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
// file          : convexity_check_2.C
// package       : Convex_hull (3.3)
// maintainer    : Stefan Schirra <stschirr@mpi-sb.mpg.de>
// source        : convex_hull_2.lw
// revision      : 3.3
// revision_date : 03 Aug 2000
// author(s)     : Stefan Schirra
//
// coordinator   : MPI, Saarbruecken
// ======================================================================


#ifndef CGAL_CONVEXITY_CHECK_2_C
#define CGAL_CONVEXITY_CHECK_2_C

#include <CGAL/Kernel/traits_aids.h>
#include <CGAL/convexity_check_2.h>
#include <CGAL/functional.h>

CGAL_BEGIN_NAMESPACE
template <class ForwardIterator, class Traits>
bool
is_ccw_strongly_convex_2( ForwardIterator first, ForwardIterator last, 
                          const Traits& ch_traits)
{
  typedef  typename Traits::Less_xy_2      Less_xy;
  typedef  typename Traits::Leftturn_2     Leftturn;

  Less_xy  smaller_xy = ch_traits.less_xy_2_object();
  Leftturn leftturn = ch_traits.leftturn_2_object();

  ForwardIterator iter1;
  ForwardIterator iter2;
  ForwardIterator iter3;

  if ( first == last) return true;

  iter2 = first;
  iter3 = ++iter2;

  if (iter3 == last ) return true;

  ++iter3;

  if (iter3 == last ) return ( *first != *iter2 );

  iter1 = first;
  short int f = 0;

  while (iter3 != last) 
  {
      if ( !leftturn( *iter1, *iter2, *iter3 ) ) return false; 
      if ( smaller_xy( *iter2, *iter1 ) && smaller_xy( *iter2, *iter3 )) ++f;

      ++iter1;
      ++iter2;
      ++iter3;
  }

  iter3 = first;
  if ( !leftturn( *iter1, *iter2, *iter3 ) ) return false; 
  if ( smaller_xy( *iter2, *iter1 ) && smaller_xy( *iter2, *iter3 )) ++f;


  iter1 = iter2;
  iter2 = first;
  ++iter3;
  if ( !leftturn( *iter1, *iter2, *iter3 ) ) return false; 
  if ( smaller_xy( *iter2, *iter1 ) && smaller_xy( *iter2, *iter3 )) ++f;


  return ( f > 1 ) ? false : true;
}

template <class ForwardIterator, class Traits>
bool
is_cw_strongly_convex_2( ForwardIterator first, ForwardIterator last, 
                         const Traits& ch_traits)
{
  typedef  typename Traits::Less_xy_2       Less_xy;
  typedef  typename Traits::Leftturn_2      Leftturn;
  typedef  Rightturn_by_leftturn< Leftturn> Rightturn;

  Less_xy  smaller_xy = ch_traits.less_xy_2_object();
  Rightturn rightturn = ch_traits.leftturn_2_object();

  ForwardIterator iter1;
  ForwardIterator iter2;
  ForwardIterator iter3;

  if ( first == last) return true;

  iter2 = first;
  iter3 = ++iter2;

  if (iter3 == last ) return true;

  ++iter3;

  if (iter3 == last ) return ( *first != *iter2 );

  iter1 = first;
  short int f = 0;

  while (iter3 != last) 
  {
      if ( !rightturn( *iter1, *iter2, *iter3 ) ) return false;
      if ( smaller_xy( *iter2, *iter1 ) && smaller_xy( *iter2, *iter3 )) ++f;

      ++iter1;
      ++iter2;
      ++iter3;
  }

  iter3 = first;
  if ( !rightturn( *iter1, *iter2, *iter3 ) ) return false;
  if ( smaller_xy( *iter2, *iter1 ) && smaller_xy( *iter2, *iter3 )) ++f;


  iter1 = iter2;
  iter2 = first;
  ++iter3;
  if ( !rightturn( *iter1, *iter2, *iter3 ) ) return false;
  if ( smaller_xy( *iter2, *iter1 ) && smaller_xy( *iter2, *iter3 )) ++f;


  return ( f > 1 ) ? false : true;
}

template <class ForwardIterator1, class ForwardIterator2, class Traits>
bool
ch_brute_force_check_2(ForwardIterator1 first1, ForwardIterator1 last1,
                            ForwardIterator2 first2, ForwardIterator2 last2,
                            const Traits&  ch_traits)
{
  typedef    typename Traits::Leftturn_2    Left_of_line;
  ForwardIterator1 iter11;
  ForwardIterator2 iter21;
  ForwardIterator2 iter22;

  if ( first1 == last1) return true;

  if ( first2 == last2) return false;

  if ( successor(first2) == last2 )
  {
      while (first1 != last1)
      {
          if ( *first1++ != *first2 ) return false;
      }
      return true;
  }

  Left_of_line  left_turn = ch_traits.leftturn_2_object();
  iter22 = first2;
  iter21 = iter22++;
  while (iter22 != last2)
  {
      iter11 = std::find_if( first1, last1, 
                             bind_1(bind_1(left_turn,*iter22++), *iter21++) );
      if (iter11 != last1 ) return false;
  }

  iter11 = std::find_if( first1, last1, 
                         bind_1(bind_1(left_turn, *first2), *iter21) );
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
  typedef    typename Traits::Left_of_line_2    Left_of_line;
  ForwardIterator1 iter11;
  ForwardIterator2 iter21;
  ForwardIterator2 iter22;

  if ( first1 == last1) return true;

  if ( first2 == last2) return false;

  if ( successor(first2) == last2 ) return true;

  Left_of_line  left_turn = ch_traits.leftturn_2_object();
  iter22 = first2;
  iter21 = iter22++;
  while (iter22 != last2)
  {
      iter11 = std::find_if( first1, last1, 
                             bind_1(bind_1(left_turn, *iter22++), *iter21++) );
      if (iter11 != last1 ) return false;
  }

  return true;
}


CGAL_END_NAMESPACE

#endif // CGAL_CONVEXITY_CHECK_2_C
