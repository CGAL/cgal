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
// package       : $CGAL_Package: $
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
#ifndef CGAL_PROTECT_IOSTREAM_H
#include <iostream.h>
#define CGAL_PROTECT_IOSTREAM_H
#endif // CGAL_PROTECT_IOSTREAM_H

#ifndef CGAL_IO_PRINT_OFF_H
#include <CGAL/IO/print_OFF.h>
#endif // CGAL_IO_PRINT_OFF_H

#ifndef CGAL_IO_SCAN_OFF_H
#include <CGAL/IO/scan_OFF.h>
#endif // CGAL_IO_SCAN_OFF_H

//  CGAL_Polyhedron
#ifndef CGAL_POLYHEDRON_3_H
#include <CGAL/Polyhedron_3.h>
#endif

template <class Traits, class HDS> inline
ostream& operator<<( ostream& out, const CGAL_Polyhedron_3<Traits,HDS>& P) {
    // writes P to `out' in PRETTY, ASCII or BINARY format
    // as the stream indicates.
    if ( CGAL_is_pretty( out))
        CGAL_print_OFF( out, P);
    else if ( CGAL_is_binary( out))
        CGAL_print_OFF( out, P, true);
    else // ASCII without comments
        CGAL_print_OFF( out, P, false, true);
    return out;
}

template <class Traits, class HDS> inline
istream& operator>>( istream& in, CGAL_Polyhedron_3<Traits,HDS>& P) {
    // reads a polyhedron from `in' and appends it to P.
    CGAL_scan_OFF( in, P);
    return in;
}
#endif // CGAL_IO_POLYHEDRON_IOSTREAM_H //
// EOF //
