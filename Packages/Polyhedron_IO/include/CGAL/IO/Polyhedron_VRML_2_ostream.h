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
// file          : include/CGAL/IO/Polyhedron_VRML_2_ostream.h
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
// Print a Polyhedron_3 in VRML 2.0 file format (.wrl)
// ============================================================================

#ifndef CGAL_IO_POLYHEDRON_VRML_2_OSTREAM_H
#define CGAL_IO_POLYHEDRON_VRML_2_OSTREAM_H 1

#include <CGAL/basic.h>
#include <CGAL/IO/VRML_2_ostream.h>
#include <CGAL/IO/File_writer_VRML_2.h>
#include <CGAL/IO/generic_print_polyhedron.h>
#include <CGAL/Polyhedron_3.h>
#include <iostream>

CGAL_BEGIN_NAMESPACE

template < class Traits,
           class Items,
#ifndef CGAL_CFG_NO_TMPL_IN_TMPL_PARAM
           template < class T, class I, class A>
#endif
           class HDS, class Alloc>
VRML_2_ostream&
operator<<( VRML_2_ostream& out,
            const Polyhedron_3<Traits,Items,HDS,Alloc>& P) {
    File_writer_VRML_2  writer;
    generic_print_polyhedron( out.os(), P, writer);
    return out;
}

CGAL_END_NAMESPACE
#endif // CGAL_IO_POLYHEDRON_VRML_2_OSTREAM_H //
// EOF //
