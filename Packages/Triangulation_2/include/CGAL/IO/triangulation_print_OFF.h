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
// release       :
// release_date  :
//
// file          : include/CGAL/IO/triangulation_print_OFF.h
// source        : web/Triangulation_2.fw
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Lutz Kettner  <kettner@inf.ethz.ch>
//
// coordinator   : Herve Bronnimann  <Herve.Bronnimann@sophia.inria.fr>
//
// ============================================================================


// Print a Triangulation<Traits> with 3d points in object file format (OFF).

#ifndef CGAL_TRIANGULATION_PRINT_OFF_H
#define CGAL_TRIANGULATION_PRINT_OFF_H 1

#ifndef CGAL_PROTECT_MAP_H
#include <map.h>
#define CGAL_PROTECT_MAP_H
#endif // CGAL_PROTECT_MAP_H

template < class Triang >
void
CGAL_triangulation_print_OFF( ostream& out, const Triang& triang,
                              bool binary = false, bool noc = false)
{
    CGAL_precondition( triang.is_valid());
    typedef typename Triang::Vertex           Vertex;
    typedef typename Triang::Vertex_iterator  Vertex_iterator;
    typedef typename Triang::Face_iterator    Face_iterator;
    // Build a map from vertex pointers to vertex indices.
    map<const Vertex*,size_t, less<const Vertex*> > mapping;
    size_t vn = 0;
    Vertex_iterator vi = triang.vertices_begin();
    for ( ; vi != triang.vertices_end(); ++vi) {
        CGAL_assertion( ! triang.is_infinite( vi));
        mapping[ &*vi] = vn;
        vn++;
    }
    CGAL_assertion( vn == triang.number_of_vertices());

    // Count finite and infinite faces.
    size_t fn  = 0;
    Face_iterator fi = triang.faces_begin();
    for ( ; fi != triang.faces_end(); ++fi) {
        CGAL_assertion( ! triang.is_infinite( fi));
        fn++;
    }
    size_t fin = triang.number_of_faces() - fn;

    CGAL_File_writer_OFF  writer( binary, noc);
    writer.header( out, vn, 3 * fn + fin, fn);

    vi = triang.vertices_begin();
    for ( ; vi != triang.vertices_end(); ++vi) {
        CGAL_assertion( ! triang.is_infinite( vi));
        writer.write_vertex(CGAL_to_double(vi->point().x()),
                            CGAL_to_double(vi->point().y()),
                            CGAL_to_double(vi->point().z()));
    }
    writer.write_facet_header();

    fi = triang.faces_begin();
    while ( fn --) {
        writer.write_facet_begin( 3);
        CGAL_assertion( mapping.find(&*(fi->vertex(0))) != mapping.end());
        CGAL_assertion( mapping.find(&*(fi->vertex(1))) != mapping.end());
        CGAL_assertion( mapping.find(&*(fi->vertex(2))) != mapping.end());
        writer.write_facet_vertex_index( mapping[ &*(fi->vertex(0))]);
        writer.write_facet_vertex_index( mapping[ &*(fi->vertex(1))]);
        writer.write_facet_vertex_index( mapping[ &*(fi->vertex(2))]);
        writer.write_facet_end();
        ++fi;
    }
    CGAL_assertion( fi == triang.faces_end());
    writer.footer();
}

#endif // CGAL_TRIANGULATION_PRINT_OFF_H //
