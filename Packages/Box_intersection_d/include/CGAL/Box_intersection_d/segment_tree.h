#ifndef CGAL_BOX_INTERSECTION_D_SEGMENT_TREE_H
#define CGAL_BOX_INTERSECTION_D_SEGMENT_TREE_H

#include <CGAL/basic.h>
#include <CGAL/Box_intersection_d/one_way_scan.h>
#include <CGAL/Box_intersection_d/modified_two_way_scan.h>

#include <algorithm>
#include <iterator>
#include <cstdlib>
#include <cmath>
#include <functional>
#include <cassert>

// original limits in g++3.2 is totally wrong. on some systems it is even missing
#include <my_limits>


CGAL_BEGIN_NAMESPACE

#define BOX_INTERSECTION_DEBUG 0

template< class RandomAccessIter, class Predicate_traits >
RandomAccessIter
median_of_three( RandomAccessIter a, RandomAccessIter b, RandomAccessIter c,
                 Predicate_traits traits, unsigned int dim )
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

template< class RandomAccessIter, class Predicate_traits >
RandomAccessIter
iterative_radon( RandomAccessIter begin, RandomAccessIter end,
                 Predicate_traits traits, unsigned int dim, int num_levels )
{
    if( num_levels < 0 )
        return begin + lrand48() % std::distance( begin, end );

    return median_of_three(
         iterative_radon( begin, end, traits, dim, num_levels - 1 ),
         iterative_radon( begin, end, traits, dim, num_levels - 1 ),
         iterative_radon( begin, end, traits, dim, num_levels - 1 ),
	 traits, dim );
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
    int levels = (int)(.91*log(((double)std::distance(begin,end))/137.0)+1);
    levels = (levels <= 0) ? 1 : levels;
    RandomAccessIter it = iterative_radon( begin, end, traits, dim, levels );
    mi = Predicate_traits::get_lo( *it, dim );
    return std::partition( begin, end,
                           typename Predicate_traits::Lo_less( mi, dim ) );
}


#if BOX_INTERSECTION_DEBUG
 static int level = -1;
 #define DUMP(msg) { \
   for( unsigned int i = level; i; --i ) \
     std::cout << "  "; \
    std::cout << msg; \
  }
#else
 #define DUMP(msg) ;
#endif


template< class ForwardIter, class Traits >
void dump_points( ForwardIter begin, ForwardIter end, Traits traits,
                  unsigned int dim ) {
    while( begin != end ) {
        std::cout << Traits::get_lo( *begin, dim ) << " ";
        ++begin;
    }
    std::cout << std::endl;
}

template< class ForwardIter, class Traits >
void dump_intervals( ForwardIter begin, ForwardIter end, Traits traits,
                     unsigned int dim ) {
    while( begin != end ) {
        std::cout << "[" << Traits::get_lo( *begin, dim ) << ","
                         << Traits::get_hi( *begin, dim ) << ") ";
        ++begin;
    }
    std::cout << std::endl;
}

template< class ForwardIter, class  Traits >
void dump_box_numbers( ForwardIter begin, ForwardIter end, Traits traits ) {
    while( begin != end ) {
        std::cout << Traits::get_id( *begin ) << " ";
        ++begin;
    }
    std::cout << std::endl;
}

template< class T >
struct Counter {
   T& value;
   Counter( T& value ) : value( value ) { ++value; }
   ~Counter() { --value; }
};

template<class RandomAccessIter,class Callback,class T,class Predicate_traits>
void segment_tree( RandomAccessIter p_begin, RandomAccessIter p_end,
                   RandomAccessIter i_begin, RandomAccessIter i_end,
                   T lo, T hi,
                   Callback& callback, Predicate_traits traits,
                   unsigned int cutoff, unsigned int dim, bool in_order )
{
    typedef typename Predicate_traits::Spanning   Spanning;
    typedef typename Predicate_traits::Lo_less    Lo_less;
    typedef typename Predicate_traits::Hi_greater Hi_greater;

    const T inf = workaround::numeric_limits< T >::inf();
    const T sup = workaround::numeric_limits< T >::sup();

#if BOX_INTERSECTION_DEBUG
    Counter<int> bla( level );
    //DUMP("----------------===========[ new node ]============-------------")
    DUMP("range: [" << lo << "," << hi << ") dim " << dim << std::endl )
    DUMP("intervals: " )
    //dump_box_numbers( i_begin, i_end, traits );
    dump_intervals( i_begin, i_end, traits, dim );
    DUMP("points: " )
    //dump_box_numbers( p_begin, p_end, traits );
    dump_points( p_begin, p_end, traits, dim );
#endif

#if SEGMENT_TREE_CHECK_INVARIANTS
    {
        // first: each point is inside segment [lo,hi)
        for( RandomAccessIter it = p_begin; it != p_end; ++it ) {
            assert( Lo_less( hi, dim )(*it) );
            assert( Lo_less( lo, dim )(*it) == false );
        }
        // second: each interval intersects segment [lo,hi)
        for( RandomAccessIter it = i_begin; it != i_end; ++it )
            assert( Hi_greater( lo, dim )(*it) && Lo_less( hi, dim )(*it) );
    }
#endif

    if( p_begin == p_end || i_begin == i_end || lo >= hi )
        return;

    if( dim == 0 )  {
        DUMP( "dim = 0. scanning ... " << std::endl )
        one_way_scan( p_begin, p_end, i_begin, i_end,
                      callback, traits, dim, in_order );
        return;
    }

    if( (unsigned int)std::distance( p_begin, p_end ) < cutoff ||
        (unsigned int)std::distance( i_begin, i_end ) < cutoff  )
    {
        DUMP( "scanning ... " << std::endl )
        modified_two_way_scan( p_begin, p_end, i_begin, i_end,
                               callback, traits, dim, in_order );
        return;
    }

    RandomAccessIter i_span_end = lo == inf || hi == sup ? i_begin :
        std::partition( i_begin, i_end, Spanning( lo, hi, dim ) );

    if( i_begin != i_span_end ) {
        DUMP( "checking spanning intervals ... " << std::endl )
        // make two calls for roots of segment tree at next level.
        segment_tree( p_begin, p_end, i_begin, i_span_end, inf, sup,
                      callback, traits, cutoff, dim - 1,  in_order );
        segment_tree( i_begin, i_span_end, p_begin, p_end, inf, sup,
                      callback, traits, cutoff, dim - 1, !in_order );
    }

    T mi;
    RandomAccessIter p_mid = split_points( p_begin, p_end, traits, dim, mi );

    if( p_mid == p_begin || p_mid == p_end )  {
        DUMP( "unable to split points! ")
        //dump_points( p_begin, p_end, traits, dim );
        DUMP( "performing modified two_way_san ... " << std::endl )
        modified_two_way_scan( p_begin, p_end, i_span_end, i_end,
                               callback, traits, dim, in_order );
        return;
    }

    RandomAccessIter i_mid;
    // separate left intervals.
    // left intervals have a low point strictly less than mi
    i_mid = std::partition( i_span_end, i_end, Lo_less( mi, dim ) );
    DUMP("->left" << std::endl )
    segment_tree( p_begin, p_mid, i_span_end, i_mid, lo, mi,
                  callback, traits, cutoff, dim, in_order );
    // separate right intervals.
    // right intervals have a high point strictly higher than mi
    i_mid = std::partition( i_span_end, i_end, Hi_greater( mi, dim ) );
    DUMP("->right"<< std::endl )
    segment_tree( p_mid, p_end, i_span_end, i_mid, mi, hi,
                  callback, traits, cutoff, dim, in_order );
}

CGAL_END_NAMESPACE

#endif
