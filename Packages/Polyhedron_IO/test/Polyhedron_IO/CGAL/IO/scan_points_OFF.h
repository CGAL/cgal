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
// package       : $CGAL_Package: $
// source        : polyhedron_io.fw
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Lutz Kettner  <kettner@inf.ethz.ch>
//
// coordinator   : Herve Bronnimann  <Herve.Bronnimann@sophia.inria.fr>
//
// scans only the points from object file format (OFF) file
// ============================================================================

#ifndef CGAL_IO_SCAN_POINTS_OFF_H
#define CGAL_IO_SCAN_POINTS_OFF_H 1
#ifndef CGAL_BASIC_H
#include <CGAL/basic.h>
#endif

#ifndef CGAL_PROTECT_STDDEF_H
#include <stddef.h>
#define CGAL_PROTECT_STDDEF_H
#endif // CGAL_PROTECT_STDDEF_H

#ifndef CGAL_IO_FILE_HEADER_OFF_H
#include <CGAL/IO/File_header_OFF.h>
#endif // CGAL_IO_FILE_HEADER_OFF_H

#ifndef CGAL_IO_FILE_SCANNER_OFF_H
#include <CGAL/IO/File_scanner_OFF.h>
#endif // CGAL_IO_FILE_SCANNER_OFF_H

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
    typedef typename Point::RT RT;
    int  i;
    for ( i = 0; i < scanner.size_of_vertices(); i++) {
        CGAL_file_scan_vertex( scanner, *points);
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
CGAL_scan_points_OFF( istream& in, Point*& points, bool verbose = false) {
    // scans the points from a polyhedral surface in OFF from `in'
    // and stores them in a dynamically allocated array in `points'.
    // The user is responsible to `delete[]' the array when appropriate.
    // Returns the past-the-end iterator for the array `points'.
    CGAL_File_scanner_OFF scanner( in, verbose);
    if ( ! in) {
        if ( verbose) {
            cerr << " " << endl;
            cerr << "CGAL_scan_points_OFF(): "
                    "input error: file format is not in OFF." << endl;
        }
        points = 0;
        return 0;
    }
    // Allocate points.
    points = new Point[ scanner.size_of_vertices()];
    return CGAL_scan_points_OFF( scanner, points);
}
#endif // CGAL_IO_SCAN_POINTS_OFF_H //
// EOF //
