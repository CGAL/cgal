// Copyright (c) 2004  Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Lutz Kettner  <kettner@mpi-sb.mpg.de>
//                 Andreas Meyer <ameyer@mpi-sb.mpg.de>

#ifndef CGAL_BOX_INTERSECTION_D_SEGMENT_TREE_H
#define CGAL_BOX_INTERSECTION_D_SEGMENT_TREE_H

#include <CGAL/license/Box_intersection_d.h>

#include <CGAL/basic.h>
#include <CGAL/Box_intersection_d/box_limits.h>

#include <boost/random/linear_congruential.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/variate_generator.hpp>

#include <algorithm>
#include <iterator>
#include <functional>
#include <cstdlib>
#include <cmath>
#include <climits>
#include <cstddef>

namespace CGAL {

namespace Box_intersection_d {

#define CGAL_BOX_INTERSECTION_DEBUG 0


template< class ForwardIter1, class ForwardIter2,
          class Callback, class Traits >
void all_pairs( ForwardIter1 p_begin, ForwardIter1 p_end,
                ForwardIter2 i_begin, ForwardIter2 i_end,
                Callback callback, Traits , bool complete_case = false)
{
    const int last_dim = Traits::dimension() - 1;
    for( ForwardIter1 p = p_begin; p != p_end; ++p ) {
        for( ForwardIter2 i = i_begin; i != i_end; ++i ) {
            if ((complete_case && Traits::id(*p) >= Traits::id(*i))
                || Traits::id(*p) == Traits::id(*i))
                continue;
            for( int dim = 0; dim <= last_dim; ++dim )
                if( !Traits::does_intersect( *p, *i, dim ) )
                    goto no_intersection1;
            callback( *p, *i );
        no_intersection1:
            ;
        }
    }

}

template< class ForwardIter, class Callback, class Traits >
void all_pairs( ForwardIter p_begin, ForwardIter p_end,
                Callback callback, Traits)
{
    const int last_dim = Traits::dimension() - 1;
    // loops actually only up to p_end-1, but we stay with the forward iterator
    // requirement and have one unnecessary but harmless additional iteration
    for( ForwardIter p = p_begin; p != p_end; ++p ) {
        ForwardIter i = p;
        ++i;
        for( ; i != p_end; ++i ) {
            for( int dim = 0; dim <= last_dim; ++dim )
                if( !Traits::does_intersect( *p, *i, dim ) )
                    goto no_intersection1;
            callback( *p, *i );
        no_intersection1:
            ;
        }
    }
}


template< class RandomAccessIter1, class RandomAccessIter2,
          class Callback, class Traits >
void one_way_scan( RandomAccessIter1 p_begin, RandomAccessIter1 p_end,
                   RandomAccessIter2 i_begin, RandomAccessIter2 i_end,
                   Callback callback, Traits, int last_dim,
                   bool in_order = true )
{
    typedef typename Traits::Compare Compare;

    // Putting a parallel sort here slows down the overall parallel algorithm
    std::sort( p_begin, p_end, Compare( 0 ) );
    std::sort( i_begin, i_end, Compare( 0 ) );

    // for each box viewed as interval i
    for( RandomAccessIter2 i = i_begin; i != i_end; ++i ) {
        // look for the first box b with i.min <= p.min
        for( ; p_begin != p_end && Traits::is_lo_less_lo( *p_begin, *i, 0 );
             ++p_begin ) {}

        // look for all boxes with p.min < i.max
        for( RandomAccessIter1 p = p_begin;
             p != p_end && Traits::is_lo_less_hi( *p, *i, 0 );
             ++p )
        {
            if( Traits::id( *p ) == Traits::id( *i ) )
                continue;
            for( int dim = 1; dim <= last_dim; ++dim )
                if( !Traits::does_intersect( *p, *i, dim ) )
                    goto no_intersection;
            if( in_order )
                callback( *p, *i );
            else
                callback( *i, *p );
        no_intersection:
            ;
        }
    }

}

template< class RandomAccessIter1, class RandomAccessIter2,
          class Callback, class Traits >
void modified_two_way_scan(
    RandomAccessIter1 p_begin, RandomAccessIter1 p_end,
    RandomAccessIter2 i_begin, RandomAccessIter2 i_end,
    Callback callback, Traits, int last_dim,
    bool in_order = true )
{
    typedef typename Traits::Compare Compare;

    std::sort( p_begin, p_end, Compare( 0 ) );
    std::sort( i_begin, i_end, Compare( 0 ) );

    // for each box viewed as interval
    while( i_begin != i_end && p_begin != p_end ) {
        if( Traits::is_lo_less_lo( *i_begin, *p_begin, 0 ) ) {
            for( RandomAccessIter1 p = p_begin;
                 p != p_end && Traits::is_lo_less_hi( *p, *i_begin, 0 );
                 ++p )
            {
                if( Traits::id( *p ) == Traits::id( *i_begin ) )
                    continue;

                for( int dim = 1; dim <= last_dim; ++dim )
                    if( !Traits::does_intersect( *p, *i_begin, dim ) )
                        goto no_intersection1;
                if( Traits::contains_lo_point( *i_begin, *p, last_dim ) ) {
                    if( in_order )
                        callback( *p, *i_begin );
                    else
                        callback( *i_begin, *p );
                }
            no_intersection1:
                ;
            }
            ++i_begin;
        } else {
            for( RandomAccessIter2 i = i_begin;
                 i != i_end && Traits::is_lo_less_hi( *i, *p_begin, 0 );
                 ++i )
            {
                if( Traits::id( *p_begin ) == Traits::id( *i ) )
                    continue;
                for( int dim = 1; dim <= last_dim; ++dim )
                    if( !Traits::does_intersect( *p_begin, *i, dim ) )
                        goto no_intersection2;
                if( Traits::contains_lo_point( *i, *p_begin, last_dim ) ) {
                    if( in_order )
                        callback( *p_begin, *i );
                    else
                        callback( *i, *p_begin );
                }
            no_intersection2:
                ;
            }
            ++p_begin;
        }
    }

}


template< class RandomAccessIter, class Predicate_traits >
RandomAccessIter
median_of_three( RandomAccessIter a, RandomAccessIter b, RandomAccessIter c,
                 Predicate_traits, int dim )
{

    if( Predicate_traits::is_lo_less_lo( *a, *b, dim ) )
        if( Predicate_traits::is_lo_less_lo( *b, *c, dim ) )
            return b;
        else if( Predicate_traits::is_lo_less_lo( *a, *c, dim ) )
            return c;
        else
            return a;
    else if( Predicate_traits::is_lo_less_lo( *a, *c, dim ) )
        return a;
    else if( Predicate_traits::is_lo_less_lo( *b, *c, dim ) )
        return c;
    else
        return b;
}


template< class RandomAccessIter, class Predicate_traits, class Generator>
class Iterative_radon {

  RandomAccessIter begin;
  Predicate_traits traits;
  int dim;
  Generator& generator;

public:

  Iterative_radon( const RandomAccessIter& begin_, const Predicate_traits& traits_,
                   int dim_, Generator& generator_)
  : begin(begin_), traits(traits_), dim(dim_), generator(generator_)
  {}

  RandomAccessIter
  operator()(int num_levels)
  {
    if( num_levels < 0 ) {
      const std::ptrdiff_t d = generator();
        return begin + d;
    }

    return median_of_three((*this)(num_levels - 1 ),
                           (*this)(num_levels - 1 ),
                           (*this)(num_levels - 1 ),
                           traits, dim );
  }
};

template< class RandomAccessIter, class Predicate_traits >
RandomAccessIter
iterative_radon( RandomAccessIter begin, RandomAccessIter end,
                 Predicate_traits traits, int dim, int num_levels )
{
  typedef typename boost::variate_generator<boost::rand48&, boost::uniform_int<std::ptrdiff_t> > Generator;
  boost::rand48 rng;
  Generator generator(rng, boost::uniform_int<std::ptrdiff_t>(0, (end-begin)-1));
  Iterative_radon<RandomAccessIter, Predicate_traits, Generator> IR(begin, traits, dim, generator);
  return IR(num_levels);
}

// returns iterator for first element in [begin,end) which does not satisfy
// the Split_Points_Predicate: [begin,mid) contains only points strictly less
// than mi. so, elements from [mid,end) are equal or higher than mi.
template< class RandomAccessIter, class Predicate_traits, class T >
RandomAccessIter
split_points( RandomAccessIter begin, RandomAccessIter end,
              Predicate_traits traits, int dim, T& mi )
{
    // magic formula
    int levels = (int)(.91*std::log(((double)std::distance(begin,end))/137.0)+1);
    levels = (levels <= 0) ? 1 : levels;
    RandomAccessIter it = iterative_radon( begin, end, traits, dim, levels );
    mi = Predicate_traits::min_coord( *it, dim );
    return std::partition( begin, end,
                           typename Predicate_traits::Lo_less( mi, dim ) );
}


#if CGAL_BOX_INTERSECTION_DEBUG
 #define CGAL_BOX_INTERSECTION_DUMP(msg) { \
   for( unsigned int i = level; i; --i ) \
     std::cout << "  "; \
    std::cout << msg; \
  }
#else
 #define CGAL_BOX_INTERSECTION_DUMP(msg) ;
#endif


template< class ForwardIter, class Traits >
void dump_points( ForwardIter begin, ForwardIter end, Traits /* traits */,
                  int dim ) {
    while( begin != end ) {
        std::cout << Traits::min_coord( *begin, dim ) << " ";
        ++begin;
    }
    std::cout << std::endl;
}

template< class ForwardIter, class Traits >
void dump_intervals( ForwardIter begin, ForwardIter end, Traits /* traits */,
                     int dim ) {
    while( begin != end ) {
        std::cout << "[" << Traits::min_coord( *begin, dim ) << ","
                         << Traits::max_coord( *begin, dim ) << ") ";
        ++begin;
    }
    std::cout << std::endl;
}

template< class ForwardIter, class  Traits >
void dump_box_numbers( ForwardIter begin, ForwardIter end, Traits /* traits */ ) {
    while( begin != end ) {
        std::cout << Traits::id( *begin ) << " ";
        ++begin;
    }
    std::cout << std::endl;
}


template<class T>
class Has_member_report
{
private:
  template<class U, U>
  class check {};

  template<class C>
  static auto f(int) -> decltype(std::declval<C>().report(0) == true, char());

  template<class C>
  static int f(...);
public:
  static const bool value = (sizeof(f<T>(0)) == sizeof(char));
};

template <class T>
CGAL_CPP17_INLINE constexpr bool Has_member_report_v = Has_member_report<T>::value;

template <typename Callback>
inline
std::enable_if_t<Has_member_report<Callback>::value, bool>
report_impl(Callback callback, int dim)
{
  return callback.report(dim);
}

template <typename Callback>
inline
std::enable_if_t<!Has_member_report<Callback>::value, bool>
report_impl(const Callback&, int)
{
  return false;
}

template <typename Callback>
inline
std::enable_if_t<Has_member_report<Callback>::value, void>
progress_impl(Callback callback, double d)
{
  callback.progress(d);
}

template <typename Callback>
inline
std::enable_if_t<!Has_member_report<Callback>::value, void>
progress_impl(const Callback&, double)
{}




template< class T >
struct Counter {
   T& value;
   Counter( T& value ) : value( value ) { ++value; }
   ~Counter() { --value; }
};

template< class RandomAccessIter1, class RandomAccessIter2,
          class Callback, class T, class Predicate_traits >
void segment_tree( RandomAccessIter1 p_begin, RandomAccessIter1 p_end,
                   RandomAccessIter2 i_begin, RandomAccessIter2 i_end,
                   T lo, T hi,
                   Callback callback, Predicate_traits traits,
                   std::ptrdiff_t cutoff, int dim, bool in_order)
{
    typedef typename Predicate_traits::Spanning   Spanning;
    typedef typename Predicate_traits::Lo_less    Lo_less;
    typedef typename Predicate_traits::Hi_greater Hi_greater;

    const T inf = box_limits< T >::inf();
    const T sup = box_limits< T >::sup();


    CGAL_STATIC_THREAD_LOCAL_VARIABLE(int, level, -1);
    Counter<int> bla( level );

    #if CGAL_BOX_INTERSECTION_DEBUG
    CGAL_BOX_INTERSECTION_DUMP("range: [" << lo << "," << hi << ") dim "
                                          << dim << std::endl )
    CGAL_BOX_INTERSECTION_DUMP("intervals: " )
    //dump_box_numbers( i_begin, i_end, traits );
    dump_intervals( i_begin, i_end, traits, dim );
    CGAL_BOX_INTERSECTION_DUMP("points: " )
    //dump_box_numbers( p_begin, p_end, traits );
    dump_points( p_begin, p_end, traits, dim );
#endif

#if CGAL_SEGMENT_TREE_CHECK_INVARIANTS
    {
        // first: each point is inside segment [lo,hi)
        for( RandomAccessIter1 it = p_begin; it != p_end; ++it ) {
            CGAL_assertion( Lo_less( hi, dim )(*it) );
            CGAL_assertion( Lo_less( lo, dim )(*it) == false );
        }
        // second: each interval intersects segment [lo,hi)
        for( RandomAccessIter2 it = i_begin; it != i_end; ++it )
            CGAL_assertion( Hi_greater(lo,dim)(*it) && Lo_less(hi,dim)(*it));
    }
#endif

    if( p_begin == p_end || i_begin == i_end || lo >= hi ){
      if(report_impl(callback, dim)){
          progress_impl(callback, 1.0 / (1 << level));
        }
        return;
    }

    if( dim == 0 )  {
        CGAL_BOX_INTERSECTION_DUMP( "dim = 0. scanning ... " << std::endl )
        one_way_scan( p_begin, p_end, i_begin, i_end,
                      callback, traits, dim, in_order );
        if(report_impl(callback,dim)){
          progress_impl(callback, 1.0 / (1 << level));
        }
        return;
    }

    if( std::distance( p_begin, p_end ) < cutoff ||
        std::distance( i_begin, i_end ) < cutoff  )
    {
        CGAL_BOX_INTERSECTION_DUMP( "scanning ... " << std::endl )
        modified_two_way_scan( p_begin, p_end, i_begin, i_end,
                               callback, traits, dim, in_order );
        if(report_impl(callback,dim)){
          progress_impl(callback, 1.0 / (1 << level));
        }
        return;
    }

    RandomAccessIter2 i_span_end = lo == inf || hi == sup ? i_begin :
        std::partition( i_begin, i_end, Spanning( lo, hi, dim ) );

    if( i_begin != i_span_end ) {
        CGAL_BOX_INTERSECTION_DUMP( "checking spanning intervals ... "
                                    << std::endl )
        // make two calls for roots of segment tree at next level.
        segment_tree( p_begin, p_end, i_begin, i_span_end, inf, sup,
                      callback, traits, cutoff, dim - 1,  in_order );
        segment_tree( i_begin, i_span_end, p_begin, p_end, inf, sup,
                      callback, traits, cutoff, dim - 1, !in_order );
    }

    T mi;
    RandomAccessIter1 p_mid = split_points( p_begin, p_end, traits, dim, mi );

    if( p_mid == p_begin || p_mid == p_end )  {
        CGAL_BOX_INTERSECTION_DUMP( "unable to split points! ")
        //dump_points( p_begin, p_end, traits, dim );
        CGAL_BOX_INTERSECTION_DUMP( "performing modified two_way_san ... "
                                     << std::endl )
        modified_two_way_scan( p_begin, p_end, i_span_end, i_end,
                               callback, traits, dim, in_order );
        if(report_impl(callback,dim)){
          progress_impl(callback, 1.0 / (1 << level));
        }
        return;
    }

    RandomAccessIter2 i_mid;
    // separate left intervals.
    // left intervals have a low point strictly less than mi
    i_mid = std::partition( i_span_end, i_end, Lo_less( mi, dim ) );
    CGAL_BOX_INTERSECTION_DUMP("->left" << std::endl )
    segment_tree( p_begin, p_mid, i_span_end, i_mid, lo, mi,
                  callback, traits, cutoff, dim, in_order );
    // separate right intervals.
    // right intervals have a high point strictly higher than mi
    i_mid = std::partition( i_span_end, i_end, Hi_greater( mi, dim ) );
    CGAL_BOX_INTERSECTION_DUMP("->right"<< std::endl )
    segment_tree( p_mid, p_end, i_span_end, i_mid, mi, hi,
                  callback, traits, cutoff, dim, in_order );
}

#if CGAL_BOX_INTERSECTION_DEBUG
 #undef CGAL_BOX_INTERSECTION_DUMP
#endif
#undef CGAL_BOX_INTERSECTION_DEBUG

} // end namespace Box_intersection_d



} //namespace CGAL

#endif
