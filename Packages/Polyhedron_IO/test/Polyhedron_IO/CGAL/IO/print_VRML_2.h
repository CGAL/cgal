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
// file          : print_VRML_2.h
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

#ifndef CGAL_IO_PRINT_VRML_2_H
#define CGAL_IO_PRINT_VRML_2_H 1
#ifndef CGAL_IO_POLYHEDRON_VRML_2_OSTREAM_H
#include <CGAL/IO/Polyhedron_VRML_2_ostream.h>
#endif // CGAL_IO_POLYHEDRON_VRML_2_OSTREAM_H

template <class Traits, class HDS>
void
CGAL_print_VRML_2( ostream& out,
                  const CGAL_Polyhedron_3<Traits,HDS>& P) {
    CGAL_VRML_2_ostream os( out);
    os << P;
}
#endif // CGAL_IO_PRINT_VRML_2_H //
// EOF //
