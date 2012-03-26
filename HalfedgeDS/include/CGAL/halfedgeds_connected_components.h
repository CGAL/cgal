// Copyright (c) 1997  
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved. 
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
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
// Author(s)     : Lutz Kettner  <kettner@mpi-sb.mpg.de>

#ifndef CGAL_HALFEDGEDS_CONNECTED_COMPONENTS_H
#define CGAL_HALFEDGEDS_CONNECTED_COMPONENTS_H 1

#include <CGAL/Unique_hash_map.h>
#include <CGAL/iterator.h>
#include <vector>

namespace CGAL {

template <class HDS, class Output_iterator, 
          class Halfedge_iterator, class Halfedge_handle>
std::size_t 
halfedgeds_connected_components( HDS hds, Output_iterator result) {
    Unique_hash_map< Halfedge_iterator, bool> hedge_map( false);
    std::size_t count = 0;
    std::vector< Halfedge_handle> hstack;
    hstack.reserve( hds.size_of_halfedges() + 3/ 4);
    Halfedge_iterator scan = hds.halfedges_begin();
    while ( scan != hds.halfedges_end()) {
        // first trace the component, then report it
        hstack.clear();
        hstack.push_back( scan);
        while ( ! hstack.empty()) {
            Halfedge_handle h = hstack.back();
            hstack.pop_back();
            if ( ! hedge_map[ h] ) {
                hedge_map[ h] = true;
                hstack.push_back( h->next());
                hstack.push_back( h->opposite());
            }
        }
        ++count;
        *result++ = scan;
        while( hedge_map[ scan])
            ++scan;
    }
    return count;
}

template <class HDS, class Output_iterator>
std::size_t 
halfedgeds_connected_components( HDS& hds, Output_iterator result) {
    typedef typename HDS::Halfedge_iterator Halfedge_iterator; 
    typedef typename HDS::Halfedge_handle   Halfedge_handle; 
    return halfedgeds_connected_components<HDS&, Output_iterator, 
        Halfedge_iterator, Halfedge_handle>( hds, result);
}

template <class HDS, class Output_iterator>
std::size_t 
halfedgeds_connected_components( const HDS& hds, Output_iterator result) {
    typedef typename HDS::Halfedge_const_iterator Halfedge_iterator; 
    typedef typename HDS::Halfedge_const_handle   Halfedge_handle; 
    return halfedgeds_connected_components<const HDS&, Output_iterator, 
        Halfedge_iterator, Halfedge_handle>( hds, result);
}

template <class HDS>
std::size_t 
halfedgeds_connected_components( HDS& hds) {
    return halfedgeds_connected_components( hds, Emptyset_iterator());
}


} //namespace CGAL
#endif // CGAL_HALFEDGEDS_CONNECTED_COMPONENTS_H //
// EOF //
