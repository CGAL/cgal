#line 12 "cgal_header.any"
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
// file          : print_VRML.h
// source        : polyhedron_io.fw
#line 26 "cgal_header.any"
                                   
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Lutz Kettner  <kettner@inf.ethz.ch>
//
// coordinator   : MPI Saarbruecken (Stefan Schirra <stschirr@mpi-sb.mpg.de>)
//
// Print a Polyhedron<Traits> in VRML file format (.wrl)
// ============================================================================
#line 303 "polyhedron_io.any"

#ifndef CGAL_PRINT_VRML_H
#define CGAL_PRINT_VRML_H 1
#line 3313 "polyhedron_io.any"
#ifndef CGAL_BASIC_H
#include <CGAL/basic.h>
#endif

#ifndef CGAL_FILE_WRITER_INVENTOR_H
#include <CGAL/IO/File_writer_inventor.h>
#endif

#ifndef CGAL_GENERIC_PRINT_POLYHEDRON_H
#include <CGAL/IO/generic_print_polyhedron.h>
#endif

#ifndef CGAL_POLYHEDRON_H
#include <CGAL/Polyhedron.h>
#endif

// Forward declarations.
class ostream;

template <class Traits>
void
CGAL_print_VRML( ostream& out,
                const CGAL_Polyhedron<Traits>& P,
                int version = 1) {
    CGAL_assertion( version == 1 || version == 2);
    CGAL_File_writer_inventor  writer( version);
    CGAL_generic_print_polyhedron( out, P, writer);
}
#line 306 "polyhedron_io.any"
 
#endif // CGAL_PRINT_VRML_H //
// EOF //
