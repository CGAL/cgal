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
// package       : $CGAL_Package: $
// source        : polyhedron_io.fw
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Lutz Kettner  <kettner@inf.ethz.ch>
//
// coordinator   : Herve Bronnimann  <Herve.Bronnimann@sophia.inria.fr>
//
// Writer for polyhedral surfaces in OpenInventor format (.iv)
// ============================================================================

#ifndef CGAL_IO_FILE_WRITER_INVENTOR_H
#define CGAL_IO_FILE_WRITER_INVENTOR_H 1
#ifndef CGAL_BASIC_H
#include <CGAL/basic.h>
#endif // CGAL_BASIC_H
#ifndef CGAL_PROTECT_IOSTREAM_H
#include <iostream.h>
#define CGAL_PROTECT_IOSTREAM_H
#endif // CGAL_PROTECT_IOSTREAM_H

class CGAL_File_writer_inventor {
    ostream*      out;
    size_t        _facets;
public:
    CGAL_File_writer_inventor() {}
    void header( ostream& o,
                 size_t vertices,
                 size_t halfedges,
                 size_t facets);
    void footer() const;

    void write_vertex( const double& x, const double& y, const double& z) {
        *out << "            " << x << ' ' << y << ' ' << z << ',' << '\n';
    }
    void write_facet_header() const;

    void write_facet_begin( size_t no) {
        *out << "            ";
    }
    void write_facet_vertex_index( size_t index) {
        *out << index << ',';
    }
    void write_facet_end() {
        *out << "-1,\n";
    }
};
#endif // CGAL_IO_FILE_WRITER_INVENTOR_H //
// EOF //
