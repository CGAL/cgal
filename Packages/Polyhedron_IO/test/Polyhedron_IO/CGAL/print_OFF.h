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
// file          : print_OFF.h
// source        : polyhedron_io.fw
#line 26 "cgal_header.fw"
                                   
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Lutz Kettner  <kettner@inf.ethz.ch>
//
// coordinator   : MPI Saarbruecken (Stefan Schirra <stschirr@mpi-sb.mpg.de>)
//
// Print a Polyhedron<Traits> in object file format (OFF)
// ============================================================================
#line 219 "polyhedron_io.fw"

#ifndef CGAL_PRINT_OFF_H
#define CGAL_PRINT_OFF_H 1
#line 613 "polyhedron_io.fw"
#ifndef CGAL_BASIC_H
#include <CGAL/basic.h>
#endif

#ifndef CGAL_FILE_WRITER_OFF_H
#include <CGAL/File_writer_OFF.h>
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
CGAL_print_OFF( ostream& out, const CGAL_Polyhedron<Traits>& P,
               bool binary = false) {
    // writes P to `out' in ASCII format or in binary format
    // if `binary == true'.
    CGAL_File_writer_OFF  writer( binary);
    CGAL_generic_print_polyhedron( out, P, writer);
}

template <class Traits>
void
CGAL_print_OFF( ostream& out, const CGAL_Polyhedron<Traits>& P,
               const CGAL_File_info& info,
               bool binary = false) {
    // writes P to `out' in ASCII format or in binary format
    // if `binary == true'. Writes additional file info for CGAL.
    CGAL_File_writer_OFF  writer( info, binary);
    CGAL_generic_print_polyhedron( out, P, writer);
}
#line 222 "polyhedron_io.fw"
 
#endif // CGAL_PRINT_OFF_H //
// EOF //
