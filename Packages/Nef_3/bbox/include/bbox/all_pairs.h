#ifndef CGAL_BBOX_ALL_PAIRS_H
#define CGAL_BBOX_ALL_PAIRS_H

#include <algorithm>


template< class RandomAccessIter, class Callback, class Traits >
void all_pairs( RandomAccessIter p_begin, RandomAccessIter p_end,
                RandomAccessIter i_begin, RandomAccessIter i_end,
                Callback& callback, Traits traits, int last_dim, bool in_order = true )
{
    for( RandomAccessIter p = p_begin; p != p_end; ++p ) {
        for( RandomAccessIter i = i_begin; i != i_end; ++i ) {
            bool does_intersect = true;
            for( unsigned int dim = 0; dim <= last_dim; ++dim )
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
