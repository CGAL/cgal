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
// revision      : $Id$
// revision_date : $Date$
// author(s)     : Lutz Kettner  <kettner@inf.ethz.ch>
//
// coordinator   : INRIA, Sophia Antipolis
//
// test stream support for CGAL
// ============================================================================

#include <CGAL/Cartesian.h>
#include <cstddef>
#include <cassert>
#include <sstream>
#include <CGAL/IO/Ostream_iterator.h>
#include <CGAL/IO/Istream_iterator.h>

typedef CGAL::Cartesian<double>                  Rep;
typedef Rep::Point_2                             Point;

int main()
{
    typedef CGAL::Ostream_iterator<Point,std::ostringstream>  IteratorO;
    typedef CGAL::Istream_iterator<Point,std::istringstream>  IteratorI;
    {
        std::ostringstream  out;
        CGAL::IO::set_ascii_mode( out);
        assert( CGAL::IO::is_ascii( out));
        out << Point( 1, 2) << '\0';
        std::istringstream in( out.str() );
        CGAL::IO::set_ascii_mode(in);
        assert( CGAL::IO::is_ascii(in));
        Point p;
        in >> p;
        assert( p == Point( 1, 2));
    }
    {
        std::ostringstream  out;
        CGAL::IO::set_ascii_mode( out);
        IteratorO   o(out);
        *o = Point( 1, 2);
        out << '\0';
        std::istringstream  in( out.str() );
        CGAL::IO::set_ascii_mode( in);
        IteratorI   i(in);
        Point p = *i;
        assert( p == Point( 1, 2));
    }
    return 0;
}
