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
// file          : print_inventor.h
// source        : polyhedron_io.fw
#line 26 "cgal_header.fw"
                                   
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Lutz Kettner  <kettner@inf.ethz.ch>
//
// coordinator   : MPI Saarbruecken (Stefan Schirra <stschirr@mpi-sb.mpg.de>)
//
// Print a Polyhedron<Traits> in Inventor file format (.iv)
// ============================================================================
#line 255 "polyhedron_io.fw"

#ifndef CGAL_PRINT_INVENTOR_H
#define CGAL_PRINT_INVENTOR_H 1
#line 2565 "polyhedron_io.fw"
#ifndef CGAL_FILE_WRITER_INVENTOR_H
#include <CGAL/File_writer_inventor.h>
#endif

#ifndef CGAL_GENERIC_PRINT_POLYHEDRON_H
#include <CGAL/generic_print_polyhedron.h>
#endif

#ifndef CGAL_POLYHEDRON_H
#include <CGAL/Polyhedron.h>
#endif

// Forward declarations.
class ostream;

template <class Traits>
void
CGAL_print_inventor( ostream& out, const CGAL_Polyhedron<Traits>& P) {
    CGAL_File_writer_inventor  writer( 0);
    CGAL_generic_print_polyhedron( out, P, writer);
}
#line 258 "polyhedron_io.fw"
 
#endif // CGAL_PRINT_INVENTOR_H //
// EOF //
