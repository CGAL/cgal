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
// file          : print_wavefront.h
// package       : $CGAL_Package: $
// source        : polyhedron_io.fw
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Lutz Kettner  <kettner@inf.ethz.ch>
//
// coordinator   : Herve Bronnimann  <Herve.Bronnimann@sophia.inria.fr>
//
// Print a Polyhedron_3 in Wavefront file format (.obj)
// ============================================================================

#ifndef CGAL_IO_PRINT_WAVEFRONT_H
#define CGAL_IO_PRINT_WAVEFRONT_H 1
#ifndef CGAL_IO_FILE_WRITER_WAVEFRONT_H
#include <CGAL/IO/File_writer_wavefront.h>
#endif // CGAL_IO_FILE_WRITER_WAVEFRONT_H

#ifndef CGAL_IO_GENERIC_PRINT_POLYHEDRON_H
#include <CGAL/IO/generic_print_polyhedron.h>
#endif // CGAL_IO_GENERIC_PRINT_POLYHEDRON_H

#ifndef CGAL_POLYHEDRON_3_H
#include <CGAL/Polyhedron_3.h>
#endif

// Forward declarations.
class ostream;

template <class Traits, class HDS>
void
CGAL_print_wavefront( ostream& out,
                     const CGAL_Polyhedron_3<Traits,HDS>& P) {
    CGAL_File_writer_wavefront  writer;
    CGAL_generic_print_polyhedron( out, P, writer);
}
#endif // CGAL_IO_PRINT_WAVEFRONT_H //
// EOF //
