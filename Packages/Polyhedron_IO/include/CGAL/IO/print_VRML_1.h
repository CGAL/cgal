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
// file          : include/CGAL/IO/print_VRML_1.h
// package       : Polyhedron_IO 2.11 (04 Feb 2000)
// chapter       : Support Library
//
// revision      : $Revision$
// revision_date : $Date$
//
// author(s)     : Lutz Kettner  <kettner@mpi-sb.mpg.de>
// maintainer    :
// coordinator   : INRIA, Sophia Antipolis
//
// Print a Polyhedron_3 in VRML 1.0 file format (.wrl)
// ============================================================================

#ifndef CGAL_IO_PRINT_VRML_1_H
#define CGAL_IO_PRINT_VRML_1_H 1

#include <CGAL/IO/Polyhedron_VRML_1_ostream.h>

CGAL_BEGIN_NAMESPACE

template <class Polyhedron>
void print_polyhedron_VRML_1( std::ostream& out, const Polyhedron& P) {
    VRML_1_ostream os( out);
    os << P;
}

// Deprecated global functions, replaced with functions above

template < class Traits,
           class Items,
#ifndef CGAL_CFG_NO_TMPL_IN_TMPL_PARAM
           template < class T, class I, class A>
#endif
           class HDS, class Alloc>
void
print_VRML_1( std::ostream& out, 
              const Polyhedron_3<Traits,Items,HDS,Alloc>& P) {
    VRML_1_ostream os( out);
    os << P;
}

CGAL_END_NAMESPACE
#endif // CGAL_IO_PRINT_VRML_1_H //
// EOF //
