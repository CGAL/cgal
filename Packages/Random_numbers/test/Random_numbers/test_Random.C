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
// release       : $CGAL_Revision: CGAL-wip $
// release_date  : $CGAL_Date$
//
// file          : test/Random_numbers/test_Random.C
// package       : Random_numbers 2.2.3 (20 Mar 2001)
// chapter       : Random Numbers Generator
//
// source        : web/Random.aw
// revision      : 2.4
// revision_date : 2000/11/01 13:48:56
//
// author(s)     : Sven Schönherr
// maintainer    : Sven Schönherr <sven@inf.ethz.ch>
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
    CGAL::Random::State  state;
    CGAL::default_random.save_state( state);
    
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
    
    // test state functions
    {
        CGAL::default_random.restore_state( state);     // `default_random' and
        CGAL::Random rnd( state);              // `rnd' have the same state now
        assert( CGAL::default_random.get_bool()         == rnd.get_bool());
        assert( CGAL::default_random.get_int( -100,100)
                                                    == rnd.get_int( -100,100));
        assert( CGAL::default_random.get_double()       == rnd.get_double());
        assert( CGAL::default_random                    == rnd);
    
        long init = CGAL::default_random( 9999);
        CGAL::Random rnd1( init), rnd2( init);
        assert( rnd1.get_bool()         == rnd2.get_bool()        );
        assert( rnd1.get_int( -100,100) == rnd2.get_int( -100,100));
        assert( rnd1.get_double()       == rnd2.get_double()      );
        assert( rnd1                    == rnd2                   );
    }

    return( 0);
}

// ===== EOF ==================================================================
