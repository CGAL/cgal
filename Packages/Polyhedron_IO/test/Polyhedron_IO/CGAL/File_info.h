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
// file          : File_info.h
// source        : polyhedron_io.fw
#line 26 "cgal_header.fw"
                                   
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Lutz Kettner  <kettner@inf.ethz.ch>
//
// coordinator   : MPI Saarbruecken (Stefan Schirra <stschirr@mpi-sb.mpg.de>)
//
// File information for CGAL.
// ============================================================================
#line 87 "polyhedron_io.fw"

#ifndef CGAL_FILE_INFO_H
#define CGAL_FILE_INFO_H 1
#line 658 "polyhedron_io.fw"
#ifndef CGAL_BASIC_H
#include <CGAL/basic.h>
#endif
#ifndef CGAL_C3_MATRIX_H
#include <CGAL/MatrixC3.h>
#endif

struct  CGAL_File_info {
    typedef  CGAL_MatrixC3<double>  Matrix;
    typedef  CGAL_VectorC3<double>  Vector;

    Matrix   linear;
    bool     normalized_to_sphere;
    double   radius;
    bool     rounded;
    int      rounded_bits;
    bool     terrain;
    Vector   terrain_vector;

    CGAL_File_info() :
        linear                ( 1.0),
        normalized_to_sphere  ( false),
        radius                ( 0.0),
        rounded               ( false),
        rounded_bits          ( 0),
        terrain               ( false),
        terrain_vector        ( 0.0, 0.0, 0.0)
    {}
};
#line 90 "polyhedron_io.fw"
  
#endif // CGAL_FILE_INFO_H //
// EOF //
