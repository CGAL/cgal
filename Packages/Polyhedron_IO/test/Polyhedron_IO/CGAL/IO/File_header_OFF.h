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
// file          : File_header_OFF.h
// package       : $CGAL_Package: $
// source        : polyhedron_io.fw
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Lutz Kettner  <kettner@inf.ethz.ch>
//
// coordinator   : Herve Bronnimann  <Herve.Bronnimann@sophia.inria.fr>
//
// File header information of an object file format (OFF) file
// ============================================================================

#ifndef CGAL_IO_FILE_HEADER_OFF_H
#define CGAL_IO_FILE_HEADER_OFF_H 1
#ifndef CGAL_IO_FILE_INFO_H
#include <CGAL/IO/File_info.h>
#endif // CGAL_IO_FILE_INFO_H

// Forward declarations.
class istream;

// Info structure for OFF file headers
// ===================================
class CGAL_File_header_OFF {
protected:
    bool _verbose;     // Print error messages if true.
    bool _skel;        // signals SKEL format instead of OFF format.
    int  n_vertices;
    int  n_halfedges;
    int  n_facets;
    int  _offset;      // index offset for vertices
    bool _colors;      // COFF detected.
    bool _normals;     // NOFF format stores also normals at vertices.
    bool _tag4;        // 4OFF detected.
    bool _tagDim;      // nOFF detected (will not be supported).
    int  _dim;         // dimension for nOFF (will not be supported).
    bool _binary;      // OFF in binary format.
    CGAL_File_info  _file_info;
public:
    CGAL_File_header_OFF( istream& in, bool verbose = false);

    bool verbose() const            { return _verbose; }
    bool is_SKEL() const            { return _skel; }   // SKEL format.
    bool is_OFF()  const            { return ! _skel; } // OFF format.
    int  size_of_vertices()  const  { return n_vertices; }
    int  size_of_halfedges() const  { return n_halfedges; }
    int  size_of_facets()    const  { return n_facets; }
    int  index_offset() const       { return _offset; }
    bool has_colors() const         { return _colors; } // COFF detected.
    bool has_normals() const        { return _normals;} // NOFF format.
    bool is_homogeneous() const     { return _tag4; }   // 4OFF detected.
                           // nOFF detected. (will not be supported).
    bool n_dimsional() const        { return _tagDim; }
                           // dimension for nOFF (will not be supported).
    int  dimension() const          { return _dim; }
    bool is_binary() const          { return _binary; } // OFF BINARY.
    CGAL_File_info&
         file_info()                { return _file_info; }
    const CGAL_File_info&
         file_info() const          { return _file_info; }
};
#endif // CGAL_IO_FILE_HEADER_OFF_H //
// EOF //
