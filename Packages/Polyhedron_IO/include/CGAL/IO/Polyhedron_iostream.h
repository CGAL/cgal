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
// file          : Polyhedron_iostream.h
// chapter       : $CGAL_Chapter: Support Library ... $
// package       : $CGAL_Package: Polyhedron_IO 2.11 (04 Feb 2000) $
// source        : polyhedron_io.fw
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Lutz Kettner  <kettner@inf.ethz.ch>
//
// coordinator   : Herve Bronnimann  <Herve.Bronnimann@sophia.inria.fr>
//
// Stream operators for Polyhedron_3 IO in object file format (OFF)
// ============================================================================

#ifndef CGAL_IO_POLYHEDRON_IOSTREAM_H
#define CGAL_IO_POLYHEDRON_IOSTREAM_H 1

#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/print_OFF.h>
#include <CGAL/IO/scan_OFF.h>
#include <iostream>

CGAL_BEGIN_NAMESPACE

#ifdef CGAL_USE_POLYHEDRON_DESIGN_ONE
template <class Traits, class HDS>
std::ostream& 
operator<<( std::ostream& out, const Polyhedron_3<Traits,HDS>& P) {
#else // CGAL_USE_POLYHEDRON_DESIGN_ONE //
template < class Traits,
           class Items,
#ifndef CGAL_CFG_NO_TMPL_IN_TMPL_PARAM
           template < class T, class I, class A>
#endif
           class HDS, class Alloc>
std::ostream&
operator<<( std::ostream& out, const Polyhedron_3<Traits,Items,HDS,Alloc>& P) {
#endif // CGAL_USE_POLYHEDRON_DESIGN_ONE //
    // writes P to `out' in PRETTY, ASCII or BINARY format
    // as the stream indicates.
    File_header_OFF header( is_binary( out), ! is_pretty( out), false);
    CGAL::print_OFF( out, P, header);
    return out;
}

#ifdef CGAL_USE_POLYHEDRON_DESIGN_ONE
template <class Traits, class HDS> inline
std::istream& operator>>( std::istream& in, Polyhedron_3<Traits,HDS>& P) {
#else // CGAL_USE_POLYHEDRON_DESIGN_ONE //
template < class Traits,
           class Items,
#ifndef CGAL_CFG_NO_TMPL_IN_TMPL_PARAM
           template < class T, class I, class A>
#endif
           class HDS, class Alloc>
std::istream& operator>>(std::istream& in,
                         Polyhedron_3<Traits,Items,HDS,Alloc>& P) {
#endif // CGAL_USE_POLYHEDRON_DESIGN_ONE //
    // reads a polyhedron from `in' and appends it to P.
    CGAL::scan_OFF( in, P);
    return in;
}

CGAL_END_NAMESPACE
#endif // CGAL_IO_POLYHEDRON_IOSTREAM_H //
// EOF //
