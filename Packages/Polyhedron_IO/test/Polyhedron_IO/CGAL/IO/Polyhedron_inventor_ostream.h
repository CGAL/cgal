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
// file          : Polyhedron_inventor_ostream.h
// package       : $CGAL_Package: $
// source        : polyhedron_io.fw
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Lutz Kettner  <kettner@inf.ethz.ch>
//
// coordinator   : Herve Bronnimann  <Herve.Bronnimann@sophia.inria.fr>
//
// Print a Polyhedron_3 in Inventor file format (.iv)
// ============================================================================

#ifndef CGAL_IO_POLYHEDRON_INVENTOR_OSTREAM_H
#define CGAL_IO_POLYHEDRON_INVENTOR_OSTREAM_H 1
#ifndef CGAL_IO_INVENTOR_OSTREAM_H
#include <CGAL/IO/Inventor_ostream.h>
#endif // CGAL_IO_INVENTOR_OSTREAM_H
#ifndef CGAL_IO_FILE_WRITER_INVENTOR_H
#include <CGAL/IO/File_writer_inventor.h>
#endif // CGAL_IO_FILE_WRITER_INVENTOR_H
#ifndef CGAL_IO_GENERIC_PRINT_POLYHEDRON_H
#include <CGAL/IO/generic_print_polyhedron.h>
#endif // CGAL_IO_GENERIC_PRINT_POLYHEDRON_H
#ifndef CGAL_POLYHEDRON_3_H
#include <CGAL/Polyhedron_3.h>
#endif

template <class Traits, class HDS>
CGAL_Inventor_ostream_base&
operator<< ( CGAL_Inventor_ostream_base& out,
             const CGAL_Polyhedron_3<Traits,HDS>& P) {
    CGAL_File_writer_inventor  writer;
    CGAL_generic_print_polyhedron( out.os(), P, writer);
    return out;
}
#endif // CGAL_IO_POLYHEDRON_INVENTOR_OSTREAM_H //
// EOF //
