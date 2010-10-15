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
// revision      : $Id$
// revision_date : $Date$
//
// author(s)     : Sven Schönherr <sven@inf.ethz.ch>
// coordinator   : INRIA Sophia-Antipolis
//
// implementation: test program for Random Numbers Generator
// ============================================================================

#include <CGAL/Random.h>
#include <cassert>
#include <iterator>

int
main()
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

    // test get_bits
    {
        int p1[2] = {0,};
        int p2[4] = {0,};
        int p3[8] = {0,};
        int p4[16] = {0,};
        for (int loops=0; loops < (1<<16); ++loops) {
          unsigned int l1 = CGAL::default_random.get_bits<1>();
          unsigned int l2 = CGAL::default_random.get_bits<2>();
          unsigned int l3 = CGAL::default_random.get_bits<3>();
          unsigned int l4 = CGAL::default_random.get_bits<4>();
	  assert(l1 < 2);
	  assert(l2 < 4);
	  assert(l3 < 8);
	  assert(l4 < 16);
	  // std::cout << l1 << " " << l2 << " "
	  //           << l3 << " " << l4 << std::endl;
	  ++(p1[l1]);
	  ++(p2[l2]);
	  ++(p3[l3]);
	  ++(p4[l4]);
	}
	std::cout << "Pseudo random distribution of get_bits<>():" << std::endl;
	std::copy(p1, p1+2, std::ostream_iterator<int>(std::cout, " "));
	std::cout << std::endl;
	std::copy(p2, p2+4, std::ostream_iterator<int>(std::cout, " "));
	std::cout << std::endl;
	std::copy(p3, p3+8, std::ostream_iterator<int>(std::cout, " "));
	std::cout << std::endl;
	std::copy(p4, p4+16, std::ostream_iterator<int>(std::cout, " "));
	std::cout << std::endl;
    }

    // test operator()
    {
        int  i = CGAL::default_random( 5555);
        assert( ( 0 <= i) && ( i < 5555));
    }

    // test get_seed()
    {
      CGAL::Random r (53);
      assert (r.get_seed() == 53);
    }

    // test save/restore state
    {
      CGAL::Random r1(17);
      CGAL::Random r2(23);
      CGAL::Random::State s;
      r1.save_state(s);
      r2.restore_state(s);
      assert (r1 == r2);
    }

    return( 0);
}

// ===== EOF ==================================================================
