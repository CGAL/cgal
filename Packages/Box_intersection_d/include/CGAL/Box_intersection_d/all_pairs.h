#ifndef CGAL_BOX_INTERSECTION_D_ALL_PAIRS_H
#define CGAL_BOX_INTERSECTION_D_ALL_PAIRS_H

#include <CGAL/basic.h>
#include <algorithm>

CGAL_BEGIN_NAMESPACE

template< class RandomAccessIter, class Callback, class Traits >
void all_pairs( RandomAccessIter p_begin, RandomAccessIter p_end,
                RandomAccessIter i_begin, RandomAccessIter i_end,
                Callback& callback, Traits traits, unsigned int last_dim )
{
    for( RandomAccessIter p = p_begin; p != p_end; ++p ) {
        for( RandomAccessIter i = i_begin; i != i_end; ++i ) {
            if (Traits::get_id(*p) >= Traits::get_id(*i) )
                continue;
            for( unsigned int dim = 0; dim <= last_dim; ++dim )
                if( !Traits::does_intersect( *p, *i, dim ) )
                    goto no_intersection1;
            callback( *p, *i );
        no_intersection1:
            ;
        }
    }

}

CGAL_END_NAMESPACE

#endif
