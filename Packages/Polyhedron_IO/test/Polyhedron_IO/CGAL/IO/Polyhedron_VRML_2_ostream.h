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
// file          : Polyhedron_VRML_2_ostream.h
// package       : $CGAL_Package: $
// source        : polyhedron_io.fw
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Lutz Kettner  <kettner@inf.ethz.ch>
//
// coordinator   : Herve Bronnimann  <Herve.Bronnimann@sophia.inria.fr>
//
// Print a Polyhedron_3 in VRML 2.0 file format (.wrl)
// ============================================================================

#ifndef CGAL_IO_POLYHEDRON_VRML_2_OSTREAM_H
#define CGAL_IO_POLYHEDRON_VRML_2_OSTREAM_H 1
#ifndef CGAL_BASIC_H
#include <CGAL/basic.h>
#endif
#ifndef CGAL_PROTECT_IOSTREAM_H
#include <iostream.h>
#define CGAL_PROTECT_IOSTREAM_H
#endif // CGAL_PROTECT_IOSTREAM_H
#ifndef CGAL_IO_VRML_2_OSTREAM_H
#include <CGAL/IO/VRML_2_ostream.h>
#endif // CGAL_IO_VRML_2_OSTREAM_H
#ifndef CGAL_IO_FILE_WRITER_VRML_2_H
#include <CGAL/IO/File_writer_VRML_2.h>
#endif // CGAL_IO_FILE_WRITER_VRML_2_H
#ifndef CGAL_IO_GENERIC_PRINT_POLYHEDRON_H
#include <CGAL/IO/generic_print_polyhedron.h>
#endif // CGAL_IO_GENERIC_PRINT_POLYHEDRON_H
#ifndef CGAL_POLYHEDRON_3_H
#include <CGAL/Polyhedron_3.h>
#endif

template <class Traits, class HDS>
CGAL_VRML_2_ostream&
operator<<( CGAL_VRML_2_ostream& out,
            const CGAL_Polyhedron_3<Traits,HDS>& P) {
    CGAL_File_writer_VRML_2  writer;
    CGAL_generic_print_polyhedron( out.os(), P, writer);
    return out;
}
#endif // CGAL_IO_POLYHEDRON_VRML_2_OSTREAM_H //
// EOF //
