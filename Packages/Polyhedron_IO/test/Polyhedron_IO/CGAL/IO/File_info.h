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
// package       : $CGAL_Package: $
// source        : polyhedron_io.fw
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Lutz Kettner  <kettner@inf.ethz.ch>
//
// coordinator   : Herve Bronnimann  <Herve.Bronnimann@sophia.inria.fr>
//
// File information for CGAL.
// ============================================================================

#ifndef CGAL_IO_FILE_INFO_H
#define CGAL_IO_FILE_INFO_H 1
#ifndef CGAL_BASIC_H
#include <CGAL/basic.h>
#endif
#ifdef CGAL_IO_FILE_INFO_H
#ifndef CGAL_CARTESIAN_H
#include <CGAL/Cartesian.h>
#endif
#ifndef CGAL_VECTOR_3_H
#include <CGAL/Vector_3.h>
#endif
#ifndef CGAL_AFF_TRANSFORMATION_3_H
#include <CGAL/Aff_transformation_3.h>
#endif
#else // CGAL_IO_FILE_INFO_H //
// This is code from the CEBAP project and useless in CGAL.
#ifndef CGAL__MATRIX_3_H
#include <CGAL/_Matrix_3.h>
#endif
#endif // CGAL_IO_FILE_INFO_H //

struct  CGAL_File_info {
#ifdef CGAL_IO_FILE_INFO_H
    typedef  CGAL_Aff_transformation_3< CGAL_Cartesian<double> >  Matrix;
    typedef  CGAL_Vector_3< CGAL_Cartesian<double> >              Vector;
#else // CGAL_IO_FILE_INFO_H //
    // This is code from the CEBAP project now useless in CGAL.
    typedef  CGAL__Matrix_3<double>  Matrix;
    typedef  CGAL_Vector_3<double>  Vector;
#endif

    Matrix   linear;
    bool     normalized_to_sphere;
    double   radius;
    bool     rounded;
    int      rounded_bits;
    bool     terrain;
    Vector   terrain_vector;

    CGAL_File_info() :
#ifdef CGAL_IO_FILE_INFO_H
        linear                ( CGAL_SCALING, 1.0),
#else // CGAL_IO_FILE_INFO_H //
        // This is code from the CEBAP project now useless in CGAL.
        linear                ( 1.0),
#endif
        normalized_to_sphere  ( false),
        radius                ( 0.0),
        rounded               ( false),
        rounded_bits          ( 0),
        terrain               ( false),
        terrain_vector        ( 0.0, 0.0, 0.0)
    {}
};
#endif // CGAL_IO_FILE_INFO_H //
// EOF //
