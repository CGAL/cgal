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
// package       : $CGAL_Package: $
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
#ifndef CGAL_IO_POLYHEDRON_SCAN_OFF_H
#include <CGAL/IO/Polyhedron_scan_OFF.h>
#endif // CGAL_IO_POLYHEDRON_SCAN_OFF_H

//  CGAL_Polyhedron
#ifndef CGAL_POLYHEDRON_3_H
#include <CGAL/Polyhedron_3.h>
#endif

template <class Traits, class HDS> inline
void CGAL_scan_OFF( istream& in, CGAL_Polyhedron_3<Traits,HDS>& P,
                   bool verbose = false) {
    // reads a polyhedron from `in' and appends it to P.
    typedef CGAL_Polyhedron_scan_OFF<HDS> Scanner;
    Scanner scanner( in, verbose);
    P.delegate(scanner);
}

template <class Traits, class HDS> inline
void CGAL_scan_OFF( istream& in,
                   CGAL_Polyhedron_3<Traits,HDS>& P,
                   CGAL_File_info& i, bool verbose = false) {
    // reads a polyhedron from `in' and appends it to P.
    // Returns also the CGAL_File_info structure of the polyhedron.
    typedef CGAL_Polyhedron_scan_OFF<HDS> Scanner;
    Scanner scanner( in, verbose);
    P.delegate(scanner);
    i = scanner.file_info();
}
#endif // CGAL_IO_SCAN_OFF_H //
// EOF //
