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
// file          : File_writer_wavefront.h
// source        : polyhedron_io.fw
#line 26 "cgal_header.fw"
                                   
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Lutz Kettner  <kettner@inf.ethz.ch>
//
// coordinator   : MPI Saarbruecken (Stefan Schirra <stschirr@mpi-sb.mpg.de>)
//
// Writer for polyhedral surfaces in Wavefront file format (.obj)
// ============================================================================
#line 165 "polyhedron_io.fw"

#ifndef CGAL_FILE_WRITER_WAVEFRONT_H
#define CGAL_FILE_WRITER_WAVEFRONT_H 1
#line 2358 "polyhedron_io.fw"
#ifndef CGAL_BINARY_FILE_IO_H
#include <CGAL/binary_file_io.h>
#endif

// Forward declarations.
class ostream;

class CGAL_File_writer_wavefront {
    ostream*  out;
    size_t   _facets;
public:
    void header( ostream& out,
                 size_t vertices,
                 size_t halfedges,
                 size_t facets);
    void footer() const {
        *out << "\n# End of Wavefront obj format #" << endl;
    }
    void write_vertex( const double& x, const double& y, const double& z) {
        *out << "v " << x << ' ' << y << ' ' << z << '\n';
    }
    void write_facet_header() {
        *out << "\n# " << _facets << " facets\n";
        *out << "# ------------------------------------------\n\n";
    }
    void write_facet_begin( size_t) {  // no unused
        *out << "f ";
    }
    void write_facet_vertex_index( size_t index) {
        *out << ' ' << index + 1;
    }
    void write_facet_end() {
        *out << '\n';
    }
};
#line 168 "polyhedron_io.fw"
  
#endif // CGAL_FILE_WRITER_WAVEFRONT_H //
// EOF //
