#ifndef CGAL_BOX_INTERSECTION_D_MODIFIED_TWO_WAY_SCAN_H
#define CGAL_BOX_INTERSECTION_D_MODIFIED_TWO_WAY_SCAN_H

#include <algorithm>

template< class RandomAccessIter,
          class Callback,
          class Traits >
void modified_two_way_scan( RandomAccessIter p_begin, RandomAccessIter p_end,
                   RandomAccessIter i_begin, RandomAccessIter i_end,
                   Callback& callback, Traits traits, unsigned int last_dim,
                   bool in_order = true )
{
    typedef typename Traits::Compare Compare;

    std::sort( p_begin, p_end, Compare( 0 ) );
    std::sort( i_begin, i_end, Compare( 0 ) );

    // for each box viewed as interval i
    while( i_begin != i_end && p_begin != p_end ) {
        if( Traits::is_lo_less_lo( *i_begin, *p_begin, 0 ) ) {
            for( RandomAccessIter p = p_begin;
                 p != p_end && Traits::is_lo_less_hi( *p, *i_begin, 0 );
                 ++p )
            {
                for( unsigned int dim = 1; dim <= last_dim; ++dim )
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
            for( RandomAccessIter i = i_begin;
                 i != i_end && Traits::is_lo_less_hi( *i, *p_begin, 0 );
                 ++i )
            {
                for( unsigned int dim = 1; dim <= last_dim; ++dim )
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


#endif
