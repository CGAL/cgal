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
// file          : CGAL_scan_points_OFF.h
// source        : polyhedron_io.fw
#line 26 "cgal_header.fw"
                                   
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Lutz Kettner  <kettner@inf.ethz.ch>
//
// coordinator   : MPI Saarbruecken (Stefan Schirra <stschirr@mpi-sb.mpg.de>)
//
// scans only the points from object file format (OFF) file
// ============================================================================
#line 153 "polyhedron_io.fw"

#ifndef CGAL_SCAN_POINTS_OFF_H
#define CGAL_SCAN_POINTS_OFF_H 1
#line 1671 "polyhedron_io.fw"
#ifndef _STDDEF_H
#include <stddef.h>
#endif

#ifndef CGAL_BASIC_H
#include <CGAL/basic.h>
#endif

#ifndef CGAL_FILE_HEADER_OFF_H
#include <CGAL/File_header_OFF.h>
#endif

#ifndef CGAL_FILE_SCANNER_OFF_H
#include <CGAL/File_scanner_OFF.h>
#endif

// Forward declarations.
class istream;

template <class OutputIterator, class Point>
OutputIterator
CGAL_scan_points_OFF( CGAL_File_scanner_OFF& scanner,
                     OutputIterator points, const Point*)
{
    // scans the points from a polyhedral surface in OFF from `scanner',
    // the OFF header has already been scanned.
    // read in all vertices
    typedef Point::RT RT;
    RT  x,  y,  z;  // Point coordinates.
    int  i;
    for ( i = 0; i < scanner.n_vertices; i++) {
        scanner.scan_vertex( x, y, z);
        *points = Point( x, y, z);
        ++ points;
        scanner.skip_to_next_vertex( i);
    }
    return points;
}
template <class OutputIterator>
OutputIterator
CGAL_scan_points_OFF( CGAL_File_scanner_OFF& scanner,
                     OutputIterator points)
{
    return CGAL_scan_points_OFF( scanner, points, value_type(points));
}

template <class Point>
Point*
CGAL_scan_points_OFF( istream& in, Point*& points) {
    // scans the points from a polyhedral surface in OFF from `in'
    // and stores them in a dynamnically allocated array in `points'.
    // The user is responsible to `delete[]' the array when appropriate.
    // Returns the past-the-end iterator for the array `points'.
    CGAL_File_scanner_OFF scanner( in);
    if ( scanner.error) {
        cerr << " " << endl;
        cerr << "CGAL_scan_points_OFF(): "
                "input error: file format is not in OFF." << endl;
        abort();
    }
    // Allocate points.
    points = new Point[ scanner.n_vertices];
    return CGAL_scan_points_OFF( scanner, points);
}
#line 156 "polyhedron_io.fw"
 
#endif // CGAL_SCAN_POINTS_OFF_H //
// EOF //
