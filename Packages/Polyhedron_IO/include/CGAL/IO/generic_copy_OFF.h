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
// file          : generic_copy_OFF.h
// chapter       : $CGAL_Chapter: Support Library ... $
// package       : $CGAL_Package: Polyhedron_IO 2.11 (04 Feb 2000) $
// source        : polyhedron_io.fw
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Lutz Kettner  <kettner@inf.ethz.ch>
//
// coordinator   : Herve Bronnimann  <Herve.Bronnimann@sophia.inria.fr>
//
// Generic copy of an object file format (OFF) file
// ============================================================================

#ifndef CGAL_IO_GENERIC_COPY_OFF_H
#define CGAL_IO_GENERIC_COPY_OFF_H 1
#ifndef CGAL_BASIC_H
#include <CGAL/basic.h>
#endif
#ifndef CGAL_PROTECT_CSTDDEF
#include <cstddef>
#define CGAL_PROTECT_CSTDDEF
#endif
#ifndef CGAL_IO_FILE_HEADER_OFF_H
#include <CGAL/IO/File_header_OFF.h>
#endif // CGAL_IO_FILE_HEADER_OFF_H
#ifndef CGAL_IO_FILE_SCANNER_OFF_H
#include <CGAL/IO/File_scanner_OFF.h>
#endif // CGAL_IO_FILE_SCANNER_OFF_H
#ifndef CGAL_PROTECT_IOSTREAM
#include <iostream>
#define CGAL_PROTECT_IOSTREAM
#endif

CGAL_BEGIN_NAMESPACE

template <class Writer>
void
generic_copy_OFF( File_scanner_OFF& scanner,
                  std::ostream& out,
                  Writer& writer) {
    std::istream& in = scanner.in();
    // scans a polyhedral surface in OFF from `in' and writes it
    // to `out' in the format provided by `writer'.
    if ( ! in) {
        if ( scanner.verbose()) {
            std::cerr << " " << std::endl;
            std::cerr << "generic_copy_OFF(): "
                         "input error: file format is not in OFF."
                      << std::endl;
        }
        return;
    }

    // Print header. Number of halfedges is only trusted if it is
    // a polyhedral surface.
    writer.write_header( out,
                         scanner.size_of_vertices(),
                         scanner.polyhedral_surface() ?
                             scanner.size_of_halfedges() : 0,
                         scanner.size_of_facets());

    // read in all vertices
    double  x,  y,  z;  // Point coordinates.
    int  i;
    for ( i = 0; i < scanner.size_of_vertices(); i++) {
        scanner.scan_vertex( x, y, z);
        writer.write_vertex( x, y, z);
        scanner.skip_to_next_vertex( i);
    }

    // read in all facets
    writer.write_facet_header();
    for ( i = 0; i < scanner.size_of_facets(); i++) {
        if ( ! in)
            return;
        Integer32 no;
        scanner.scan_facet( no, i);
        writer.write_facet_begin( no);
        for ( int j = 0; j < no; j++) {
            Integer32 index;
            scanner.scan_facet_vertex_index( index, i);
            writer.write_facet_vertex_index( index);
        }
        writer.write_facet_end();
        scanner.skip_to_next_facet( i);
    }
    writer.write_footer();
}

template <class Writer>
void
generic_copy_OFF( std::istream& in, std::ostream& out, Writer& writer,
                  bool verbose = false) {
    // scans a polyhedral surface in OFF from `in' and writes it
    // to `out' in the format provided by `writer'.
    File_scanner_OFF scanner( in, verbose);
    generic_copy_OFF( scanner, out, writer);
}

CGAL_END_NAMESPACE
#endif // CGAL_IO_GENERIC_COPY_OFF_H //
// EOF //
