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
// file          : Polyhedron_scan_OFF.h
// source        : polyhedron_io.fw
#line 26 "cgal_header.fw"
                                   
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Lutz Kettner  <kettner@inf.ethz.ch>
//
// coordinator   : MPI Saarbruecken (Stefan Schirra <stschirr@mpi-sb.mpg.de>)
//
// Scanner for OFF as polyhedral surface builder
// ============================================================================
#line 231 "polyhedron_io.fw"

#ifndef CGAL_POLYHEDRON_SCAN_OFF_H
#define CGAL_POLYHEDRON_SCAN_OFF_H 1
#line 1409 "polyhedron_io.fw"
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

//  CGAL_Polyhedron_builder
//  CGAL_Polyhedron_incremental_builder
#ifndef CGAL_POLYHEDRON_H
#include <CGAL/Polyhedron.h>
#endif


// Forward declarations.
class istream;

template < class R>
class CGAL_Polyhedron_scan_OFF :  public CGAL_Polyhedron_builder<R> {
protected:
    istream&       _in;
    CGAL_File_info  _file_info;
public:

// DEFINITION
//
// CGAL_Polyhedron_scan_OFF<Traits> is a polyhedral surface builder, see
// Section .... It scans a polyhedron given in OFF from a stream and
// appends it incrementally using the incremental builder from Section ...
//
// `Traits' is the representation class of the polyhedral surface to be
// constructed.
//
// `#include <CGAL/polyhedron_io.h>'
//
// TYPES

// New creation variable is: `scanner'
//
// CREATION

    CGAL_Polyhedron_scan_OFF( istream& in) : _in(in) {}
        // creates the scanner and stores the `istream' in its internal
        // state.

// Activation

    void build( Traits& poly);
        // scans the `istream' known from creation time and builds the
        // polyhedral surface. If the stream does not contain a valid OFF
        // object or does this not describe a permissable polyhedral
        // surface (e.g. a non-manifold) the scanner reports an approprate
        // error message to `cerr' and terminates the program.
        // Postcondition: the `poly' is a valid polyhedral surface.

    CGAL_File_info  file_info() const   { return _file_info; }
};

template < class Traits >
void
CGAL_Polyhedron_scan_OFF<Traits>:: build( Traits& target) {
    CGAL_Polyhedron_incremental_builder<Traits> B( target);
    CGAL_File_scanner_OFF scanner( _in);
    if ( scanner.error) {
        cerr << " " << endl;
        cerr << "CGAL_Polyhedron_scan_OFF<Traits>::" << endl;
        cerr << "build(): input error: file format is not in OFF." << endl;
        abort();
    }
    _file_info = scanner.file_info;

    B.begin_surface( scanner.n_vertices,
                     scanner.n_facets,
                     scanner.n_halfedges);

    typedef Traits::Point Point;
    typedef Point::RT  RT;
    RT  x,  y,  z;  // Point coordinates.

    // read in all vertices
    int  i;
    for ( i = 0; i < scanner.n_vertices; i++) {
        scanner.scan_vertex( x, y, z);
        B.add_vertex( Point( x, y, z));
        scanner.skip_to_next_vertex( i);
    }

    // read in all facets
    for ( i = 0; i < scanner.n_facets; i++) {
        B.begin_facet();
        CGAL_Integer32 no;
        scanner.scan_facet( no, i);
        for ( int j = 0; j < no; j++) {
            CGAL_Integer32 index;
            scanner.scan_facet_vertex_index( index, i);
            B.add_vertex_to_facet( index);
        }
        B.end_facet();
        scanner.skip_to_next_facet( i);
    }
#if 0
    B.check_unconnected_vertices(true);
    if ( ! B.remove_unconnected_vertices()) {
        cerr << " " << endl;
        cerr << "CGAL_Polyhedron_scan_OFF<Traits>::" << endl;
        cerr << "build(): input error: the unconnected vertices detected "
                "in the file cannot be removed." << endl;
        abort();
    }
#endif
    B.end_surface();
}
#line 234 "polyhedron_io.fw"
 
#endif // CGAL_POLYHEDRON_SCAN_OFF_H //
// EOF //
