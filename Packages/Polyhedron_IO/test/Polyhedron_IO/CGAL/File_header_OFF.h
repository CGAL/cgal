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
// file          : File_header_OFF.h
// source        : polyhedron_io.fw
#line 26 "cgal_header.fw"
                                   
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Lutz Kettner  <kettner@inf.ethz.ch>
//
// coordinator   : MPI Saarbruecken (Stefan Schirra <stschirr@mpi-sb.mpg.de>)
//
// File header information of an object file format (OFF) file
// ============================================================================
#line 99 "polyhedron_io.fw"

#ifndef CGAL_FILE_HEADER_OFF_H
#define CGAL_FILE_HEADER_OFF_H 1
#line 690 "polyhedron_io.fw"
#ifndef CGAL_FILE_INFO_H
#include <CGAL/File_info.h>
#endif

// Forward declarations.
class istream;

// Info structure for OFF file headers
// ===================================
struct CGAL_File_header_OFF {
    bool error;       // signals wrong file format
    bool skel;        // signals SKEL format instead of OFF format.
    int  n_vertices;
    int  n_halfedges;
    int  n_facets;
    int  offset;      // index offset for vertices
    bool colors;      // COFF detected.
    bool normals;     // NOFF format stores also normals at vertices.
    bool tag4;        // 4OFF detected.
    bool tagDim;      // nOFF detected (will not be supported).
    int  dim;         // dimension for nOFF (will not be supported).
    bool binary;      // OFF in binary format.
    CGAL_File_info  file_info;
    CGAL_File_header_OFF( istream& in);
};
#line 102 "polyhedron_io.fw"
  
#endif // CGAL_FILE_HEADER_OFF_H //
// EOF //
