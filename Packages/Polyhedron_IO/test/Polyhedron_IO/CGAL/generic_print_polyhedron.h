#line 12 "cgal_header.fw"
// ============================================================================
//
// Copyright (c) 1997 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       : $CGAL_Revision: $
// release_date  : $CGAL_Date: $
//
// file          : generic_print_polyhedron.h
// source        : polyhedron_io.fw
#line 26 "cgal_header.fw"
                                   
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Lutz Kettner  <kettner@inf.ethz.ch>
//
// coordinator   : MPI Saarbruecken (Stefan Schirra <stschirr@mpi-sb.mpg.de>)
//
// A generic writer for polyhedral surfaces parameterized by file format
// ============================================================================
#line 54 "polyhedron_io.fw"

#ifndef CGAL_GENERIC_PRINT_POLYHEDRON_H
#define CGAL_GENERIC_PRINT_POLYHEDRON_H 1
#line 413 "polyhedron_io.fw"
#ifndef CGAL_BASIC_H
#include <CGAL/basic.h>
#endif

#ifndef CGAL_POLYHEDRON_H
#include <CGAL/Polyhedron.h>
#endif

// Forward declarations.
class ostream;

template <class Traits, class Writer>
void
CGAL_generic_print_polyhedron( ostream& out,
                              const CGAL_Polyhedron<Traits>& P,
                              Writer& writer) {
    // writes P to `out' in the format provided by `writer'.
    typedef CGAL_Polyhedron<Traits>                    Poly;
    typedef Poly::Vertex_const_iterator            VCI;
    typedef Poly::Facet_const_iterator             FCI;
    typedef Poly::Halfedge_facet_const_circulator  HFCC;
    // Print header.
    writer.header( out,
                   P.size_of_vertices(),
                   P.size_of_halfedges(),
                   P.size_of_facets());
    for( VCI vi = P.vertices_begin(); vi != P.vertices_end(); ++vi) {
        writer.write_vertex( (*vi).point().x(),
                             (*vi).point().y(),
                             (*vi).point().z());
    }
    typedef CGAL_Inverse_index< VCI, Poly::Vertex, Poly::Size> Index;
    Index index( P.vertices_begin(), P.vertices_end());
    writer.write_facet_header();

    for( FCI fi = P.facets_begin(); fi != P.facets_end(); ++fi) {
        HFCC hc = (*fi).facet_begin();
        HFCC hc_end = hc;
        size_t n = CGAL_circulator_size( hc);
        CGAL_assertion( n >= 3);
        writer.write_facet_begin( n);
        do {
            writer.write_facet_vertex_index( index[ VCI((*hc).vertex())]);
            ++hc;
        } while( hc != hc_end);
        writer.write_facet_end();
    }
    writer.footer();
}
#line 57 "polyhedron_io.fw"
 
#endif // CGAL_GENERIC_PRINT_POLYHEDRON_H //
// EOF //
