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
// file          : test_support.C
// chapter       : $CGAL_Chapter: Stream Support $
// package       : $CGAL_Package: Stream_support 2.4 (28 Jul 1999) $
// source        : support.fw
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Lutz Kettner  <kettner@inf.ethz.ch>
//
// coordinator   : INRIA, Sophia Antipolis
//
// test stream support for CGAL
// ============================================================================


#ifndef CGAL_BASIC_H
#include <CGAL/basic.h>
#endif // CGAL_BASIC_H
#ifndef CGAL_CARTESIAN_H
#include <CGAL/Cartesian.h>
#endif // CGAL_CARTESIAN_H
#ifndef CGAL_PROTECT_CSTDDEF
#include <cstddef>
#define CGAL_PROTECT_CSTDDEF
#endif

#ifndef CGAL_POINT_2_H
#include <CGAL/Point_2.h>
#endif // CGAL_POINT_2_H
#ifndef CGAL_PROTECT_STRSTREAM
#include <strstream>
#define CGAL_PROTECT_STRSTREAM
#endif // CGAL_PROTECT_STRSTREAM
#ifndef CGAL_IO_OSTREAM_ITERATOR_H
#include <CGAL/IO/Ostream_iterator.h>
#endif // CGAL_IO_OSTREAM_ITERATOR_H
#ifndef CGAL_IO_ISTREAM_ITERATOR_H
#include <CGAL/IO/Istream_iterator.h>
#endif // CGAL_IO_ISTREAM_ITERATOR_H

using namespace CGAL;

void test_stream_iterator() {
    typedef Cartesian<double>                        Rep;
    typedef Point_2<Rep>                             Point;
    typedef Ostream_iterator<Point,std::ostrstream>  IteratorO;
    typedef Istream_iterator<Point,std::istrstream>  IteratorI;
    {
        char buffer[1000];
        std::ostrstream  out( buffer, 1000);
        set_ascii_mode( out);
        CGAL_assertion( is_ascii( out));
        out << Point( 1, 2) << '\0';
        std::istrstream in( buffer, 1000);
        set_ascii_mode(in);
        CGAL_assertion( is_ascii(in));
        Point p;
        in >> p;
        CGAL_assertion( p == Point( 1, 2));
    }
    {
        char buffer[1000];
        std::ostrstream  out( buffer, 1000);
        set_ascii_mode( out);
        IteratorO   o(out);
        *o = Point( 1, 2);
        out << '\0';
        std::istrstream  in( buffer, 1000);
        set_ascii_mode( in);
        IteratorI   i(in);
        Point p = *i;
        CGAL_assertion( p == Point( 1, 2));
    }
}

int main() {
    test_stream_iterator();
    return 0;
}
// EOF //
