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
// file          : File_writer_inventor.h
// source        : polyhedron_io.fw
#line 26 "cgal_header.fw"
                                   
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Lutz Kettner  <kettner@inf.ethz.ch>
//
// coordinator   : MPI Saarbruecken (Stefan Schirra <stschirr@mpi-sb.mpg.de>)
//
// Writer for polyhedral surfaces in object file format (.iv)
// ============================================================================
#line 186 "polyhedron_io.fw"

#ifndef CGAL_FILE_WRITER_INVENTOR_H
#define CGAL_FILE_WRITER_INVENTOR_H 1
#line 2455 "polyhedron_io.fw"
// Forward declarations.
class ostream;

class CGAL_File_writer_inventor {
    ostream*      out;
    size_t        _facets;
    int           _version;
    const char*   indent;
public:
    CGAL_File_writer_inventor( int version = 0) : _version( version) {
                              // version == 0: Inventor
                              // version == 1: VRML 1.0
                              // version == 2: VRML 2.0
        CGAL_assertion( version >= 0 && version <= 2);
    }
    void header( ostream& o,
                 size_t vertices,
                 size_t halfedges,
                 size_t facets);
    void footer() const;

    void write_vertex( const double& x, const double& y, const double& z) {
        *out << indent << x << ' ' << y << ' ' << z << ',' << '\n';
    }
    void write_facet_header() {
        if ( _version < 2) {
            *out << "  ]\n"
                    "}\n"
                    "# " << _facets << " facets\n"
                    "IndexedFaceSet {\n"
                    "  coordIndex [\n";
        } else {
            *out << "      ]\n"
                    "    }\n"
                    "    # " << _facets << " facets\n"
                    "    coordIndex [\n";
            indent = "      ";
        }
    }
    void write_facet_begin( size_t no) {
        *out << indent;
    }
    void write_facet_vertex_index( size_t index) {
        *out << index << ' ';
    }
    void write_facet_end() {
        *out << "0,\n";
    }
};
#line 189 "polyhedron_io.fw"
  
#endif // CGAL_FILE_WRITER_INVENTOR_H //
// EOF //
