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
// file          : scan_OFF.h
// source        : polyhedron_io.fw
#line 26 "cgal_header.fw"
                                   
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Lutz Kettner  <kettner@inf.ethz.ch>
//
// coordinator   : MPI Saarbruecken (Stefan Schirra <stschirr@mpi-sb.mpg.de>)
//
// Scan a Polyhedron<Traits> in object file format (OFF)
// ============================================================================
#line 301 "polyhedron_io.fw"

#ifndef CGAL_SCAN_OFF_H
#define CGAL_SCAN_OFF_H 1
#line 1533 "polyhedron_io.fw"
#ifndef CGAL_POLYHEDRON_SCAN_OFF_H
#include <CGAL/Polyhedron_scan_OFF.h>
#endif

//  CGAL_Polyhedron
#ifndef CGAL_POLYHEDRON_H
#include <CGAL/Polyhedron.h>
#endif

template <class Traits> inline
void CGAL_scan_OFF( istream& in, CGAL_Polyhedron<Traits>& P) {
    // reads a polyhedron from `in' and appends it to P.
    typedef CGAL_Polyhedron_scan_OFF<Traits> Scanner;
    Scanner scanner( in);
    P.delegate(scanner);
}

template <class Traits> inline
void CGAL_scan_OFF( istream& in, CGAL_Polyhedron<Traits>& P, CGAL_File_info& i) {
    // reads a polyhedron from `in' and appends it to P.
    // Returns also the CGAL_File_info structure of the polyhedron.
    typedef CGAL_Polyhedron_scan_OFF<Traits> Scanner;
    Scanner scanner( in);
    P.delegate(scanner);
    i = scanner.file_info();
}
#line 304 "polyhedron_io.fw"
 
#endif // CGAL_SCAN_OFF_H //
// EOF //
