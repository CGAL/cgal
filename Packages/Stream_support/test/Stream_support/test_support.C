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
// release       : 
// release_date  : 
//
// file          : test_support.C
// package       : Stream_support
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Lutz Kettner  <kettner@inf.ethz.ch>
//
// coordinator   : INRIA, Sophia Antipolis
//
// test stream support for CGAL
// ============================================================================

#include <CGAL/Cartesian.h>
#include <cstddef>
#include <strstream>
#include <CGAL/IO/Ostream_iterator.h>
#include <CGAL/IO/Istream_iterator.h>

using namespace CGAL;

int main()
{
    typedef Cartesian<double>                        Rep;
    typedef Rep::Point_2                             Point;
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
    return 0;
}
