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
// file          : CGAL_generic_copy_OFF.h
// source        : polyhedron_io.fw
#line 26 "cgal_header.fw"
                                   
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Lutz Kettner  <kettner@inf.ethz.ch>
//
// coordinator   : MPI Saarbruecken (Stefan Schirra <stschirr@mpi-sb.mpg.de>)
//
// Generic copy of an object file format (OFF) file
// ============================================================================
#line 141 "polyhedron_io.fw"

#ifndef CGAL_GENERIC_COPY_OFF_H
#define CGAL_GENERIC_COPY_OFF_H 1
#line 1603 "polyhedron_io.fw"
#ifndef _STDDEF_H
#include <stddef.h>
#endif

#ifndef CGAL_BASIC_H
#include <CGAL/basic.h>
#endif

#ifndef CGAL_FILE_HEADER_OFF_H
#include <CGAL/File_header_OFF.h>
#endif

#ifndef CGAL_FILE_SCANNER_OFF_H
#include <CGAL/File_scanner_OFF.h>
#endif

// Forward declarations.
class ostream;
class istream;

template <class Writer>
void
CGAL_generic_copy_OFF( istream& in, ostream& out, Writer& writer) {
    // scans a polyhedral surface in OFF from `in' and writes it
    // to `out' in the format provided by `writer'.
    CGAL_File_scanner_OFF scanner( in);
    if ( scanner.error) {
        cerr << " " << endl;
        cerr << "CGAL_generic_copy_OFF(): "
                "input error: file format is not in OFF." << endl;
        abort();
    }

    // Print header. Number of halfedges is not trusted.
    writer.header( out, scanner.n_vertices, 0, scanner.n_facets);

    // read in all vertices
    double  x,  y,  z;  // Point coordinates.
    int  i;
    for ( i = 0; i < scanner.n_vertices; i++) {
        scanner.scan_vertex( x, y, z);
        writer.write_vertex( x, y, z);
        scanner.skip_to_next_vertex( i);
    }

    // read in all facets
    writer.write_facet_header();
    for ( i = 0; i < scanner.n_facets; i++) {
        CGAL_Integer32 no;
        scanner.scan_facet( no, i);
        writer.write_facet_begin( no);
        for ( int j = 0; j < no; j++) {
            CGAL_Integer32 index;
            scanner.scan_facet_vertex_index( index, i);
            writer.write_facet_vertex_index( index);
        }
        writer.write_facet_end();
        scanner.skip_to_next_facet( i);
    }
    writer.footer();
}
#line 144 "polyhedron_io.fw"
 
#endif // CGAL_GENERIC_COPY_OFF_H //
// EOF //
