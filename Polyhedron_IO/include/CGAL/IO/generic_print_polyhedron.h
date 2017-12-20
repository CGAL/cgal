// Copyright (c) 1997  ETH Zurich (Switzerland).
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
// Author(s)     : Lutz Kettner  <kettner@mpi-sb.mpg.de>

#ifndef CGAL_IO_GENERIC_PRINT_POLYHEDRON_H
#define CGAL_IO_GENERIC_PRINT_POLYHEDRON_H 1

#include <CGAL/license/Polyhedron.h>


#include <CGAL/basic.h>
#include <CGAL/Inverse_index.h>
#include <iostream>

namespace CGAL {

template <class Polyhedron, class Writer>
void
generic_print_polyhedron( std::ostream&     out, 
                          const Polyhedron& P,
                          Writer&           writer) {
    // writes P to `out' in the format provided by `writer'.
    typedef typename Polyhedron::Vertex_const_iterator                  VCI;
    typedef typename Polyhedron::Facet_const_iterator                   FCI;
    typedef typename Polyhedron::Halfedge_around_facet_const_circulator HFCC;
    // Print header.
    writer.write_header( out,
                         P.size_of_vertices(),
                         P.size_of_halfedges(),
                         P.size_of_facets());
    for( VCI vi = P.vertices_begin(); vi != P.vertices_end(); ++vi) {
        writer.write_vertex( ::CGAL::to_double( vi->point().x()),
                             ::CGAL::to_double( vi->point().y()),
                             ::CGAL::to_double( vi->point().z()));
    }
    typedef Inverse_index< VCI> Index;
    Index index( P.vertices_begin(), P.vertices_end());
    writer.write_facet_header();

    for( FCI fi = P.facets_begin(); fi != P.facets_end(); ++fi) {
        HFCC hc = fi->facet_begin();
        HFCC hc_end = hc;
        std::size_t n = circulator_size( hc);
        CGAL_assertion( n >= 3);
        writer.write_facet_begin( n);
        do {
            writer.write_facet_vertex_index( index[ VCI(hc->vertex())]);
            ++hc;
        } while( hc != hc_end);
        writer.write_facet_end();
    }
    writer.write_footer();
}

} //namespace CGAL
#endif // CGAL_IO_GENERIC_PRINT_POLYHEDRON_H //
// EOF //
