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
// file          : test/Random_numbers/test_Random.cpp
// package       : $CGAL_Package: Random_numbers $
// chapter       : Random Numbers Generator
//
// revision      : $Id$
// revision_date : $Date$
//
// author(s)     : Sven Schönherr <sven@inf.ethz.ch>
// coordinator   : INRIA Sophia-Antipolis
//
// implementation: test program for Random Numbers Generator
// ============================================================================

#include <CGAL/Random.h>
#include <CGAL/algorithm.h>
#include <cassert>
#include <iterator>
#include <vector>
#include <algorithm>

int
main()
{
    // test get_bool
    {
      bool b = CGAL::get_default_random().get_bool();
        assert( ! b || b);
    }

    // test get_int
    {
        int  l = CGAL::get_default_random().get_int( -100, 0);
        int  u = CGAL::get_default_random().get_int( 0, 1000);
        int  i = CGAL::get_default_random().get_int( l, u);
        assert( ( l <= i) && ( i < u));

        {
          std::size_t l = 0, u = 10;
          std::size_t i = CGAL::get_default_random().uniform_int(l,u);
          assert( ( l <= i) && ( i <= u));
        }


        {
          std::ptrdiff_t l = 0, u = 10;
          std::ptrdiff_t i = CGAL::get_default_random().uniform_int(l,u);
          assert( ( l <= i) && ( i <= u));
        }

        {
          std::ptrdiff_t l = 0, u = 10;
          std::ptrdiff_t i = CGAL::get_default_random().uniform_smallint(l,u);
          assert( ( l <= i) && ( i <= u));
        }


    }

    // test get_double
    {
        double  l = CGAL::get_default_random().get_double( -123.45, -0.99);
        double  u = CGAL::get_default_random().get_double( 22.0/7.0, 33.3);
        double  d = CGAL::get_default_random().get_double( l, u);
        assert( ( l <= d) && ( d < u));

        double ho = CGAL::get_default_random().get_double(0.5);
        assert( (0.5 <= ho) && (ho < 1.0));
        double zo = CGAL::get_default_random().get_double();
        assert( (0 <= zo) && (zo < 1.0));
    }

    // test uniform_real
    {
      double d = CGAL::get_default_random().uniform_real<double>(-10.0, 10.0);
      assert( (d >= -10.0) && (d < 10.0) );

      d = CGAL::get_default_random().uniform_real<double>(0.2);
      assert( (d >= 0.2) && (d < 1.0) );

      d = CGAL::get_default_random().uniform_real<double>();
      assert( (d >= 0) && (d < 1) );

      d = CGAL::get_default_random().uniform_01<double>();
      assert( (d >= 0) && (d < 1) );
    }
  {
      float d = CGAL::get_default_random().uniform_real<float>(-10.0f, 10.0f);
      assert( (d >= -10.0f) && (d < 10.0f) );
      d = CGAL::get_default_random().uniform_real<float>(0.2f);
      assert( (d >= 0.2) && (d < 1.0) );
      d = CGAL::get_default_random().uniform_real<float>();
      assert( (d >= 0) && (d < 1) );
      d = CGAL::get_default_random().uniform_01<float>();
      assert( (d >= 0) && (d < 1) );
  }

    // test get_bits
    {
        int p1[2] = {0,};
        int p2[4] = {0,};
        int p3[8] = {0,};
        int p4[16] = {0,};
        for (int loops=0; loops < (1<<16); ++loops) {
          unsigned int l1 = CGAL::get_default_random().get_bits<1>();
          unsigned int l2 = CGAL::get_default_random().get_bits<2>();
          unsigned int l3 = CGAL::get_default_random().get_bits<3>();
          unsigned int l4 = CGAL::get_default_random().get_bits<4>();
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
      int  i = CGAL::get_default_random()( 5555);
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

    std::vector<int> numbers;
    numbers.push_back(1);
    CGAL::cpp98::random_shuffle(numbers.begin(), numbers.end(), CGAL::get_default_random());
    return( 0);
}

// ===== EOF ==================================================================
