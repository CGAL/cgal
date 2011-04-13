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
// chapter       : $CGAL_Chapter: Support Library ... $
// package       : $CGAL_Package: Polyhedron_IO 2.11 (04 Feb 2000) $
// source        : polyhedron_io.fw
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Lutz Kettner  <kettner@inf.ethz.ch>
//
// coordinator   : Herve Bronnimann  <Herve.Bronnimann@sophia.inria.fr>
//
// Print a Polyhedron_3 in object file format (OFF)
// ============================================================================

#ifndef CGAL_IO_PRINT_OFF_H
#define CGAL_IO_PRINT_OFF_H 1

#include <CGAL/basic.h>
#include <CGAL/IO/File_writer_OFF.h>
#include <CGAL/IO/generic_print_polyhedron.h>
#include <CGAL/Polyhedron_3.h>
#include <iostream>

CGAL_BEGIN_NAMESPACE

#ifdef CGAL_USE_POLYHEDRON_DESIGN_ONE
template <class Traits, class HDS>
void
print_OFF( std::ostream& out,
           const Polyhedron_3<Traits,HDS>& P,
           bool verbose = false) {
#else // CGAL_USE_POLYHEDRON_DESIGN_ONE //
template < class Traits,
           class Items,
#ifndef CGAL_CFG_NO_TMPL_IN_TMPL_PARAM
           template < class T, class I, class A>
#endif
           class HDS, class Alloc>
void
print_OFF( std::ostream& out,
           const Polyhedron_3<Traits,Items,HDS,Alloc>& P,
           bool verbose = false) {
#endif // CGAL_USE_POLYHEDRON_DESIGN_ONE //
    File_writer_OFF  writer( verbose);
    writer.header().set_binary( is_binary( out));
    writer.header().set_no_comments( ! is_pretty( out));
    writer.header().set_polyhedral_surface( true);
    writer.header().set_halfedges( P.size_of_halfedges());
    generic_print_polyhedron( out, P, writer);
}

#ifdef CGAL_USE_POLYHEDRON_DESIGN_ONE
template <class Traits, class HDS>
void
print_OFF( std::ostream& out,
           const Polyhedron_3<Traits,HDS>& P,
           const File_header_OFF& header) {
#else // CGAL_USE_POLYHEDRON_DESIGN_ONE //
template < class Traits,
           class Items,
#ifndef CGAL_CFG_NO_TMPL_IN_TMPL_PARAM
           template < class T, class I, class A>
#endif
           class HDS, class Alloc>
void
print_OFF( std::ostream& out,
           const Polyhedron_3<Traits,Items,HDS,Alloc>& P,
           const File_header_OFF& header) {
#endif // CGAL_USE_POLYHEDRON_DESIGN_ONE //
    File_writer_OFF  writer( header);
    writer.header().set_polyhedral_surface( true);
    writer.header().set_halfedges( P.size_of_halfedges());
    generic_print_polyhedron( out, P, writer);
}

CGAL_END_NAMESPACE
#endif // CGAL_IO_PRINT_OFF_H //
// EOF //
