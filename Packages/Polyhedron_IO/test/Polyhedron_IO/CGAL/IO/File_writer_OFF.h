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
// file          : File_writer_OFF.h
// package       : $CGAL_Package: $
// source        : polyhedron_io.fw
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Lutz Kettner  <kettner@inf.ethz.ch>
//
// coordinator   : Herve Bronnimann  <Herve.Bronnimann@sophia.inria.fr>
//
// Writer for polyhedral surfaces in object file format (OFF)
// ============================================================================

#ifndef CGAL_IO_FILE_WRITER_OFF_H
#define CGAL_IO_FILE_WRITER_OFF_H 1
#ifndef CGAL_IO_BINARY_FILE_IO_H
#include <CGAL/IO/binary_file_io.h>
#endif // CGAL_IO_BINARY_FILE_IO_H
#ifndef CGAL_IO_FILE_INFO_H
#include <CGAL/IO/File_info.h>
#endif // CGAL_IO_FILE_INFO_H

// Forward declarations.
class ostream;

class CGAL_File_writer_OFF {
    ostream*  out;
    size_t   _facets;
    bool     _binary;
    bool     _nocomments;
    bool     _skel;
    const CGAL_File_info*  file_info;
public:
    CGAL_File_writer_OFF( bool binary = false, bool nocomments = false,
                         bool skel = false)
        : _binary( binary),
          _nocomments( nocomments),
          _skel( skel),
          file_info(0) {}
    CGAL_File_writer_OFF( const CGAL_File_info& info,
                         bool binary = false, bool nocomments = false)
        : _binary( binary),
          _nocomments( nocomments),
          _skel( false),
          file_info( &info) {}
    void set_skel( bool skel)   { _skel = skel; }
    void header( ostream& out,
                 size_t vertices,
                 size_t halfedges,
                 size_t facets,
                 int    normals = false);
    void footer() {
        if ( ! _binary && ! _nocomments)
            *out << "\n\n# End of OFF #";
        *out << endl;
    }
    void write_vertex( const double& x, const double& y, const double& z) {
        if ( _binary) {
            CGAL__Binary_write_float32( *out, x);
            CGAL__Binary_write_float32( *out, y);
            CGAL__Binary_write_float32( *out, z);
        } else {
            *out << '\n' << x << ' ' << y << ' ' << z;
        }
    }
    void write_normal( const double& x, const double& y, const double& z) {
        if ( _binary) {
            CGAL__Binary_write_float32( *out, x);
            CGAL__Binary_write_float32( *out, y);
            CGAL__Binary_write_float32( *out, z);
        } else {
            *out << ' ' << ' ' << x << ' ' << y << ' ' << z;
        }
    }
    void write_facet_header() {
        if ( ! _binary) {
            if ( _nocomments)
                *out << '\n';
            else {
                *out << "\n\n# " << _facets << " facets\n";
                *out << "# ------------------------------------------\n\n";
            }
        }
    }
    void write_facet_begin( size_t no) {
        if ( _binary)
            CGAL__Binary_write_integer32( *out, no);
        else
            *out << no << ' ';
    }
    void write_facet_vertex_index( size_t index) {
        if ( _binary)
            CGAL__Binary_write_integer32( *out, index);
        else
            *out << ' ' << index;
    }
    void write_facet_end() {
        if ( _binary)
            CGAL__Binary_write_integer32( *out, 0);
        else
            *out << '\n';
    }
};
#endif // CGAL_IO_FILE_WRITER_OFF_H //
// EOF //
