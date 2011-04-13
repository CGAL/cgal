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
// file          : scan_OFF.h
// chapter       : $CGAL_Chapter: Support Library ... $
// package       : $CGAL_Package: Polyhedron_IO 2.11 (04 Feb 2000) $
// source        : polyhedron_io.fw
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Lutz Kettner  <kettner@inf.ethz.ch>
//
// coordinator   : Herve Bronnimann  <Herve.Bronnimann@sophia.inria.fr>
//
// Scan a Polyhedron_3 in object file format (OFF)
// ============================================================================

#ifndef CGAL_IO_SCAN_OFF_H
#define CGAL_IO_SCAN_OFF_H 1

#include <CGAL/IO/Polyhedron_scan_OFF.h>
#include <CGAL/Polyhedron_3.h>
#include <iostream>

CGAL_BEGIN_NAMESPACE

#ifdef CGAL_USE_POLYHEDRON_DESIGN_ONE

template <class Traits, class HDS> inline
void scan_OFF( std::istream& in, Polyhedron_3<Traits,HDS>& P,
               bool verbose = false) {
    // reads a polyhedron from `in' and appends it to P.
    typedef Polyhedron_scan_OFF<HDS> Scanner;
    Scanner scanner( in, verbose);
    P.delegate(scanner);
}

template <class Traits, class HDS> inline
void scan_OFF( std::istream& in,
               Polyhedron_3<Traits,HDS>& P,
               File_header_OFF& header) {
    // reads a polyhedron from `in' and appends it to P.
    // Returns also the File_header_OFF structure of the object.
    typedef Polyhedron_scan_OFF<HDS> Scanner;
    Scanner scanner( in, header.verbose());
    P.delegate(scanner);
    header = scanner.header();
}

#else // CGAL_USE_POLYHEDRON_DESIGN_ONE //

template < class Traits,
           class Items,
#ifndef CGAL_CFG_NO_TMPL_IN_TMPL_PARAM
           template < class T, class I, class A>
#endif
           class HDS, class Alloc>
void scan_OFF( std::istream& in,
               Polyhedron_3<Traits,Items,HDS,Alloc>& P,
               File_header_OFF& header) {
    // reads a polyhedron from `in' and appends it to P.
    // Returns also the File_header_OFF structure of the object.
    typedef Polyhedron_3<Traits,Items,HDS,Alloc> Polyhedron;
    typedef typename Polyhedron::HalfedgeDS HalfedgeDS;
    typedef Polyhedron_scan_OFF<HalfedgeDS> Scanner;
    Scanner scanner( in, header.verbose());
    P.delegate(scanner);
    header = scanner.header();
}

template < class Traits,
           class Items,
#ifndef CGAL_CFG_NO_TMPL_IN_TMPL_PARAM
           template < class T, class I, class A>
#endif
           class HDS, class Alloc>
void scan_OFF( std::istream& in, Polyhedron_3<Traits,Items,HDS,Alloc>& P,
               bool verbose = false) {
    // reads a polyhedron from `in' and appends it to P.
    typedef Polyhedron_3<Traits,Items,HDS,Alloc> Polyhedron;
    typedef typename Polyhedron::HalfedgeDS HalfedgeDS;
    typedef Polyhedron_scan_OFF<HalfedgeDS> Scanner;
    Scanner scanner( in, verbose);
    P.delegate(scanner);
}

#endif // CGAL_USE_POLYHEDRON_DESIGN_ONE //


CGAL_END_NAMESPACE
#endif // CGAL_IO_SCAN_OFF_H //
// EOF //
