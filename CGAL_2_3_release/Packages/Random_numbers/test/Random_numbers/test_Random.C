// ============================================================================
//
// Copyright (c) 1997-2001 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       : $CGAL_Revision: CGAL-I $
// release_date  : $CGAL_Date$
//
// file          : test/Random_numbers/test_Random.C
// package       : $CGAL_Package: Random_numbers $
// chapter       : Random Numbers Generator
//
// source        : web/Random.aw
// revision      : $Revision$
// revision_date : $Date$
//
// author(s)     : Sven Schönherr <sven@inf.ethz.ch>
// coordinator   : INRIA Sophia-Antipolis
//
// implementation: test program for Random Numbers Generator
// ============================================================================

// includes
#include <CGAL/Random.h>
#include <cassert>

int
main( int, char**)
{
    // test get_bool
    {
        bool b = CGAL::default_random.get_bool();
        assert( ! b || b);
    }
    
    // test get_int
    {
        int  l = CGAL::default_random.get_int( -100, 0);
        int  u = CGAL::default_random.get_int( 0, 1000);
        int  i = CGAL::default_random.get_int( l, u);
        assert( ( l <= i) && ( i < u));
    }
    
    // test get_double
    {
        double  l = CGAL::default_random.get_double( -123.45, -0.99);
        double  u = CGAL::default_random.get_double( 22.0/7.0, 33.3);
        double  d = CGAL::default_random.get_double( l, u);
        assert( ( l <= d) && ( d < u));
    }
    
    // test operator()
    {
        int  i = CGAL::default_random( 5555);
        assert( ( 0 <= i) && ( i < 5555));
    }

    return( 0);
}

// ===== EOF ==================================================================
