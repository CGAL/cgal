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
// chapter       : $CGAL_Chapter: Support Library ... $
// package       : $CGAL_Package: Polyhedron_IO 2.11 (04 Feb 2000) $
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

#include <CGAL/IO/Inventor_ostream.h>
#include <CGAL/IO/File_writer_inventor.h>
#include <CGAL/IO/generic_print_polyhedron.h>
#include <CGAL/Polyhedron_3.h>

CGAL_BEGIN_NAMESPACE

#ifdef CGAL_USE_POLYHEDRON_DESIGN_ONE
template <class Traits, class HDS>
Inventor_ostream_base&
operator<<( Inventor_ostream_base& out, const Polyhedron_3<Traits,HDS>& P) {
#else // CGAL_USE_POLYHEDRON_DESIGN_ONE //
template < class Traits,
           class Items,
#ifndef CGAL_CFG_NO_TMPL_IN_TMPL_PARAM
           template < class T, class I, class A>
#endif
           class HDS, class Alloc>
Inventor_ostream_base&
operator<<( Inventor_ostream_base& out, 
            const Polyhedron_3<Traits,Items,HDS,Alloc>& P) {
#endif // CGAL_USE_POLYHEDRON_DESIGN_ONE //
    File_writer_inventor  writer;
    generic_print_polyhedron( out.os(), P, writer);
    return out;
}

CGAL_END_NAMESPACE
#endif // CGAL_IO_POLYHEDRON_INVENTOR_OSTREAM_H //
// EOF //

