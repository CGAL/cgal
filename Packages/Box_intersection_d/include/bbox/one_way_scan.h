#ifndef CGAL_BBOX_ONE_WAY_SCAN_H
#define CGAL_BBOX_ONE_WAY_SCAN_H

#include <algorithm>


template< class RandomAccessIter, class Callback, class Traits >
void one_way_scan( RandomAccessIter p_begin, RandomAccessIter p_end,
                   RandomAccessIter i_begin, RandomAccessIter i_end,
                   Callback& callback, Traits traits, unsigned int last_dim,
                   bool in_order = true )
{
    typedef typename Traits::Compare Compare;
    std::sort( p_begin, p_end, Compare( 0 ) );
    std::sort( i_begin, i_end, Compare( 0 ) );

    // for each box viewed as interval i
    for( RandomAccessIter i = i_begin; i != i_end; ++i ) {
        // look for the first box b with i.min <= p.min
        for( ; p_begin != p_end && Traits::is_lo_less_lo( *p_begin, *i, 0 ); 
             ++p_begin );

        // look for all boxes with p.min < i.max
        for( RandomAccessIter p = p_begin;
             p != p_end && Traits::is_lo_less_hi( *p, *i, 0 );
             ++p )
        {
            for( unsigned int dim = 1; dim <= last_dim; ++dim )
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


#endif
