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
// package       : $CGAL_Package: $
// source        : polyhedron_io.fw
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Lutz Kettner  <kettner@inf.ethz.ch>
//
// coordinator   : Herve Bronnimann  <Herve.Bronnimann@sophia.inria.fr>
//
// Scanner for OFF as polyhedral surface builder
// ============================================================================

#ifndef CGAL_IO_POLYHEDRON_SCAN_OFF_H
#define CGAL_IO_POLYHEDRON_SCAN_OFF_H 1
#ifndef CGAL_BASIC_H
#include <CGAL/basic.h>
#endif

#ifndef CGAL_PROTECT_STDDEF_H
#include <stddef.h>
#define CGAL_PROTECT_STDDEF_H
#endif // CGAL_PROTECT_STDDEF_H

#ifndef CGAL_IO_FILE_HEADER_OFF_H
#include <CGAL/IO/File_header_OFF.h>
#endif // CGAL_IO_FILE_HEADER_OFF_H

#ifndef CGAL_IO_FILE_SCANNER_OFF_H
#include <CGAL/IO/File_scanner_OFF.h>
#endif // CGAL_IO_FILE_SCANNER_OFF_H

#ifndef CGAL_MODIFIER_BASE_H
#include <CGAL/Modifier_base.h>
#endif

#ifndef CGAL_POLYHEDRON_INCREMENTAL_BUILDER_3_H
#include <CGAL/Polyhedron_incremental_builder_3.h>
#endif // CGAL_POLYHEDRON_INCREMENTAL_BUILDER_3_H

#ifndef CGAL_POLYHEDRON_3_H
#include <CGAL/Polyhedron_3.h>
#endif


// Forward declarations.
class istream;

template < class HDS>
class CGAL_Polyhedron_scan_OFF :  public CGAL_Modifier_base<HDS> {
protected:
    istream&       _in;
    bool           _verbose;
    CGAL_File_info  _file_info;
public:

    typedef HDS Halfedge_data_structure;

// DEFINITION
//
// CGAL_Polyhedron_scan_OFF<Traits> is a polyhedral surface builder.
// It scans a polyhedron given in OFF from a stream and appends it
// incrementally using the incremental builder.

    CGAL_Polyhedron_scan_OFF( istream& in, bool verbose = false)
        : _in(in), _verbose( verbose) {}

// Activation

    void operator()( HDS& hds);

    CGAL_File_info  file_info() const   { return _file_info; }
};

template < class HDS >
void
CGAL_Polyhedron_scan_OFF<HDS>:: operator()( HDS& target) {
    CGAL_File_scanner_OFF scanner( _in, _verbose);
    if ( ! _in) {
        if ( _verbose) {
            cerr << " " << endl;
            cerr << "CGAL_Polyhedron_scan_OFF<HDS>::" << endl;
            cerr << "operator(): input error: file format is not in OFF."
                 << endl;
        }
        return;
    }
    _file_info = scanner.file_info();

    CGAL_Polyhedron_incremental_builder_3<HDS> B( target, _verbose);
    B.begin_surface( scanner.size_of_vertices(),
                     scanner.size_of_facets(),
                     scanner.size_of_halfedges());

    typedef typename HDS::Point Point;
    typedef typename Point::RT  RT;

    // read in all vertices
    int  i;
    for ( i = 0; i < scanner.size_of_vertices(); i++) {
        Point p;
        CGAL_file_scan_vertex( scanner, p);
        B.add_vertex( p);
        scanner.skip_to_next_vertex( i);
    }
    if ( ! _in  || B.error()) {
        B.rollback();
        _in.clear( ios::badbit);
        return;
    }

    // read in all facets
    for ( i = 0; i < scanner.size_of_facets(); i++) {
        B.begin_facet();
        CGAL_Integer32 no;
        scanner.scan_facet( no, i);
        if( ! _in || B.error() || no < 3) {
            if ( _verbose) {
                cerr << " " << endl;
                cerr << "CGAL_Polyhedron_scan_OFF<Traits>::" << endl;
                cerr << "operator()(): input error: facet " << i
                     << " has less than 3 vertices." << endl;
            }
            B.rollback();
            _in.clear( ios::badbit);
            return;
        }
        for ( int j = 0; j < no; j++) {
            CGAL_Integer32 index;
            scanner.scan_facet_vertex_index( index, i);
            B.add_vertex_to_facet( index);
        }
        B.end_facet();
        scanner.skip_to_next_facet( i);
    }
    if ( ! _in  || B.error()) {
        B.rollback();
        _in.clear( ios::badbit);
        return;
    }
#if 0
    B.check_unconnected_vertices();
    if ( ! B.remove_unconnected_vertices()) {
        if ( _verbose) {
            cerr << " " << endl;
            cerr << "CGAL_Polyhedron_scan_OFF<HDS>::" << endl;
            cerr << "operator(): input error: the unconnected vertices "
                    "detected in the file cannot be removed." << endl;
        }
        B.rollback();
        _in.clear( ios::badbit);
        return;
    }
#endif
    B.end_surface();
}
#endif // CGAL_IO_POLYHEDRON_SCAN_OFF_H //
// EOF //
