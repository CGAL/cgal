// Copyright (c) 2003  INRIA Sophia-Antipolis (France) and
//                     Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// Authors : Monique Teillaud <Monique.Teillaud@sophia.inria.fr>
//           Sylvain Pion     <Sylvain.Pion@sophia.inria.fr>
// 
// Partially supported by the IST Programme of the EU as a Shared-cost
// RTD (FET Open) Project under Contract No  IST-2000-26473 
// (ECG - Effective Computational Geometry for Curves and Surfaces) 

// Test program for Root_of_2.

#include <iostream>
#include <cassert>

#include <CGAL/Random.h>

#include <CGAL/MP_Float.h>
#include <CGAL/Quotient.h>
#include <CGAL/Root_of_2.h>

#ifdef CGAL_USE_GMP
#  include <CGAL/Gmpz.h>
#  include <CGAL/Gmpq.h>
#endif

#ifdef CGAL_USE_GMPXX
#  include <CGAL/gmpxx.h>
#endif

#ifdef CGAL_USE_LEDA
#  include <CGAL/leda_real.h>
#endif

#ifdef CGAL_USE_CORE
#  include <CGAL/CORE_Expr.h>
#endif

// We should put these in a nested namespace and use it.
using CGAL::compare; // for double
using CGAL::sign;    // for double

CGAL::Random rnd;

int test_loops;

// Generate random integral values of 7 bits.
template < typename NT >
NT my_rand()
{
  return rnd.get_int(-64, 63);
}

// Generate random Root_of_2 of degree 1
template < typename Root >
Root my_rand_root_1()
{
  typedef typename Root::RT  RT;
  Root r;
  do {
    RT a = my_rand<RT>();
    if (a == 0)
      continue;
    r = Root(a, my_rand<RT>());
  } while (! is_valid(r));
  return r;
}

// Generate random Root_of_2 of degree 2
template < typename Root >
Root my_rand_root()
{
  typedef typename Root::RT  RT;
  do {
    RT a = my_rand<RT>();
    RT b = my_rand<RT>();
    RT c = my_rand<RT>();
    if (a == 0 && b == 0)
      return Root(c);
    if (a == 0)
      return Root(b, c);
    if (b*b-4*a*c >= 0)
      return Root(a, b, c, rnd.get_bool());
  } while (true);
}

// Compares the comparison results of two roots with
// the comparison results of their double approximations.
template < typename Root >
bool
test_compare(const Root &r1, const Root &r2)
{
  if (compare(r1, r2) == compare(to_double(r1), to_double(r2)))
    return true;

  std::cout << "ERROR in comparing :" << std::endl << " r1 = ";
  print(std::cout, r1);
  std::cout << " [approx = " << to_double(r1) << " ]\n" << " r2 = ";
  print(std::cout, r2);
  std::cout << " [approx = " << to_double(r2) << " ]\n"
	    << " comparison gives : " << compare(r1, r2) << " for exact, but "
	    << compare(to_double(r1), to_double(r2))
	    << " for double approximation" << std::endl;
  return false;
}

template < typename Root >
bool
test_to_interval(const Root &r1)
{
  std::pair<double, double> the_interval = to_interval(r1);
  return (to_double(r1) >= the_interval.first) && (to_double(r1) <= the_interval.second);
}



template < typename Root >
bool
test_root_of()
{
  typedef typename Root::RT RT;

  std::cout << "  Testing zeros" << std::endl;
  Root zero1(1, 0, 0, true);
  Root zero2(1, 0, 0, false);
  Root zero3(0);

  assert(zero1.discriminant() == 0);
  assert(zero2.discriminant() == 0);
  assert(zero1.conjugate() == zero2);
  assert(zero2.conjugate() == zero1);
  assert(zero3.conjugate() == zero3);
  assert(is_valid(zero1));
  assert(is_valid(zero2));
  assert(is_valid(zero3));
  assert(to_double(zero1) == 0.0);
  assert(to_double(zero2) == 0.0);
  assert(to_double(zero3) == 0.0);
  assert(test_to_interval(zero1));
  assert(test_to_interval(zero2));
  assert(test_to_interval(zero3));
  assert(sign(zero1) == 0);
  assert(sign(zero2) == 0);
  assert(sign(zero3) == 0);
  assert(compare(zero1, zero1) == 0);
  assert(compare(zero2, zero2) == 0);
  assert(compare(zero1, zero2) == 0);
  assert(compare(zero2, zero1) == 0);
  assert(zero2 == zero1);
  assert(zero2 == zero3);
  assert(zero1 == zero3);
  assert(zero2 <= zero1);
  assert(zero2 <= zero3);
  assert(zero1 <= zero3);
  assert(zero2 >= zero1);
  assert(zero2 >= zero3);
  assert(zero1 >= zero3);
  assert(!(zero2 != zero1));
  assert(!(zero2 != zero3));
  assert(!(zero1 != zero3));
  assert(!(zero2 < zero1));
  assert(!(zero2 < zero3));
  assert(!(zero1 < zero3));
  assert(!(zero2 > zero1));
  assert(!(zero2 > zero3));
  assert(!(zero1 > zero3));

  assert(compare(zero2, -zero1) == 0);
  assert(compare(zero1, -zero2) == 0);

  assert(compare(zero1, zero3) == 0);
  assert(compare(zero2, zero3) == 0);

  std::cout << "  Testing ones" << std::endl;
  Root one1 (-1, 0, 1, false);
  Root mone1(-1, 0, 1, true);
  Root one2  = Root(1);
  Root mone2 = Root(-1);
  Root one3  = Root(1, -1);
  Root mone3 = Root(1, 1);
  
  assert(one1.conjugate() == mone1);
  assert(one2.conjugate() == mone2);
  assert(one3.conjugate() == mone3);
  assert(one1.discriminant()  == 4);
  assert(mone1.discriminant() == 4);
  assert(is_valid(one1));
  assert(is_valid(mone1));
  assert(is_valid(one2));
  assert(is_valid(mone2));
  assert(to_double(one1)  ==  1.0);
  assert(to_double(mone1) == -1.0);
  assert(to_double(one2)  ==  1.0);
  assert(to_double(mone2) == -1.0);
  assert(test_to_interval(one1));
  assert(test_to_interval(mone1));
  assert(test_to_interval(one2));
  assert(test_to_interval(mone2));
  assert(test_to_interval(one3));
  assert(test_to_interval(mone3));

  assert(sign( one1) > 0);
  assert(sign( one2) > 0);
  assert(sign(mone1) < 0);
  assert(sign(mone2) < 0);
  assert(compare( one1,  one1) == 0);
  assert(compare(mone1, mone1) == 0);
  assert(compare(mone1,  one1) < 0);
  assert(compare( one1, mone1) > 0);

  assert(compare( one2,  one2) == 0);
  assert(compare(mone2, mone2) == 0);
  assert(compare(mone2,  one2) < 0);
  assert(compare( one2, mone2) > 0);

  assert(compare( one1, -mone1) == 0);
  assert(compare(-one1,  mone1) == 0);

  assert(compare( one1,  one2) == 0);
  assert(compare(mone1, mone2) == 0);
  assert(compare( one1,  one3) == 0);
  assert(compare(mone1, mone3) == 0);
  assert(compare( one2,  one3) == 0);
  assert(compare(mone2, mone3) == 0);

  assert(compare( one1, zero1) > 0);
  assert(compare(mone1, zero2) < 0);

  std::cout << "  Testing output" << std::endl;
  std::cout << "  zero = ";
  print(std::cout, zero1);
  std::cout << " approx =  " << zero1 << std::endl;

  std::cout << "  zero = ";
  print(std::cout, Root(0));
  std::cout << " approx =  " << Root(0) << std::endl;

  std::cout << "  one  = ";
  print(std::cout, one1);
  std::cout << " approx =  " << one1 << std::endl;

  std::cout << "  mone = ";
  print(std::cout, mone1);
  std::cout << " approx =  " << mone1 << std::endl;

  std::cout << "  Testing degree 1" << std::endl;
  assert(compare( Root(0),  Root(0)) == 0);
  assert(compare( Root(1),  Root(1)) == 0);
  assert(compare(Root(-1), Root(-1)) == 0);
  assert(compare( Root(0),  Root(1)) < 0);
  assert(compare( Root(1), Root(-1)) > 0);

  assert(compare( Root(2, 4), Root(-2)) == 0);
  assert(compare( Root(2, 4), Root(-1)) < 0);
  assert(compare( Root(2, 4), Root(-3)) > 0);

  std::cout << "  Testing degree 1 and 2" << std::endl;
  Root rone(1, 0, -1, false);
  Root rmone(1, 0, -1, true);
  assert(rone.conjugate() == rmone);
  assert(test_to_interval(rone));
  assert(test_to_interval(rmone));
  // Compare the two roots of the above polynomial (-1 and 1),
  // succesively with -2, -1, 0, 1, 2.
  assert(compare(Root(-2), rmone) < 0);
  assert(compare(Root(-2), rone)  < 0);
  assert(compare(rmone, Root(-2)) > 0);
  assert(compare(rone, Root(-2))  > 0);
  assert(compare(Root(-1), rmone) == 0);
  assert(compare(Root(-1), rone)  < 0);
  assert(compare(rmone, Root(-1)) == 0);
  assert(compare(rone, Root(-1))  > 0);
  assert(compare(Root(0), rmone)  > 0);
  assert(compare(Root(0), rone)   < 0);
  assert(compare(rmone, Root(0))  < 0);
  assert(compare(rone, Root(0))   > 0);
  assert(compare(Root(1), rmone)  > 0);
  assert(compare(Root(1), rone)   == 0);
  assert(compare(rmone, Root(1))  < 0);
  assert(compare(rone, Root(1))   == 0);
  assert(compare(Root(2), rmone)  > 0);
  assert(compare(Root(2), rone)   > 0);
  assert(compare(rmone, Root(2))  < 0);
  assert(compare(rone, Root(2))   < 0);



  std::cout << "  Testing random roots constructed by constants" << std::endl;
  for (int i = 0; i < test_loops; ++i) {
    Root r1 = Root(my_rand<RT>());
    Root r2 = Root(my_rand<RT>());
    assert(r1 == r1);
    assert(test_to_interval(r1));
    assert(r2 == r2);
    assert(test_to_interval(r2));
    assert(test_compare(r1, r2));
  }

  std::cout << "  Testing random roots of degree 1" << std::endl;
  for (int i = 0; i < test_loops; ++i) {
    Root r1 = my_rand_root_1<Root>();
    Root r2 = my_rand_root_1<Root>();
    assert(r1 == r1);
    assert(test_to_interval(r1));
    assert(r2 == r2);
    assert(test_to_interval(r2));
    assert(test_compare(r1, r2));
  }

  std::cout << "  Testing random roots of degree 2" << std::endl;
  for (int i = 0; i < test_loops; ++i) {
    Root r1 = my_rand_root<Root>();
    Root r2 = my_rand_root<Root>();
    assert(r1 == r1);
    assert(test_to_interval(r1));
    assert(r2 == r2);
    assert(test_to_interval(r2));
    assert(test_compare(r1, r2));
  }

  std::cout << "  Testing random roots of degree 1 and 2" << std::endl;
  for (int i = 0; i < test_loops; ++i) {
    Root r1 = my_rand_root_1<Root>();
    Root r2 = my_rand_root<Root>();
    assert(r1 == r1);
    assert(test_to_interval(r1));
    assert(r2 == r2);
    assert(test_to_interval(r2));
    assert(test_compare(r1, r2));
  }

  std::cout << "  Testing squares of random roots of degree 2" << std::endl;
  for (int i = 0; i < test_loops; ++i) {
    Root r1 = my_rand_root<Root>();
    Root r2 = my_rand_root<Root>();
    Root r1_sqr = square(r1);
    Root r2_sqr = square(r2);
    assert(test_to_interval(r1));
    assert(test_to_interval(r2));
    assert(test_to_interval(r1_sqr));
    assert(test_to_interval(r2_sqr));
    double dr1 = to_double(r1);
    double dr2 = to_double(r2);
    assert(compare(r1_sqr, r2_sqr) == compare(dr1*dr1, dr2*dr2));
  }

  std::cout << "  Testing addition/subtraction of Root_of_2<NT> with NT"
            << std::endl;
  for (int i = 0; i < test_loops; ++i) {
    RT r    = my_rand<RT>();
    Root r1 = my_rand_root<Root>();
    Root r2 = r1 - r;
    Root r3 = r1 + r;
    Root r4 = r - r1;
    Root r5 = r + r1;
    assert(r1 == r1);
    assert(r2 == r2);
    assert(r3 == r3);
    assert(r4 == r4);
    assert(r5 == r5);
    assert(test_to_interval(r1));
    assert(test_to_interval(r2));
    assert(test_to_interval(r3));
    assert(test_to_interval(r4));
    assert(test_to_interval(r5));
    assert(test_compare(r1, r2));
    assert(test_compare(r1, r3));
    assert(test_compare(r1, r4));
    assert(test_compare(r1, r5));
    assert(CGAL_NTS compare(r1, r2) ==   (int) CGAL_NTS sign(r));
    assert(CGAL_NTS compare(r1, r3) == - (int) CGAL_NTS sign(r));
    assert(CGAL_NTS compare(r1, r5) == - (int) CGAL_NTS sign(r));
    assert(CGAL_NTS compare(r, r4)  == (int) CGAL_NTS sign(r1));
  }

  std::cout << "  Testing multiplication of Root_of_2<NT> with NT" << std::endl;
  for (int i = 0; i < test_loops; ++i) {
    RT r    = my_rand<RT>();
    Root r1 = my_rand_root<Root>();
    Root r2 = r1 * r;
    Root r3 = r * r1;
    assert(r1 == r1);
    assert(r2 == r2);
    assert(r3 == r3);
    assert(test_to_interval(r1));
    assert(test_to_interval(r2));
    assert(test_to_interval(r3));
    assert(test_compare(r1, r2));
    assert(test_compare(r1, r3));
    assert(test_compare(r2, r3));
    assert(r2 == r3);
    if (r > 0) {
      assert(compare(r1, r2) == compare(1, r) * sign(r1));
      assert(compare(r1, r3) == compare(1, r) * sign(r1));
    } else if (r < 0) {
      assert(compare(r1, -r2) == compare(1, -r) * sign(r1));
      assert(compare(r1, -r3) == compare(1, -r) * sign(r1));
    } else {
      assert(r2 == Root(0));
      assert(r3 == Root(0));
    }
  }

  return true;
}

int main(int argc, char **argv) {

  using CGAL::Root_of_2;

  test_loops = argc == 2 ? atoi(argv[1]) : 1000;

  bool result = true;

//  std::cout << "Testing Root_of_2<double>" << std::endl;
//  result = result && test_root_of<Root_of_2<double> >();

  std::cout << "Testing Root_of_2<MP_Float>" << std::endl;
  result = result && test_root_of<Root_of_2<CGAL::MP_Float> >();

  std::cout << "Testing Root_of_2<Quotient<MP_Float> >" << std::endl;
  result = result &&
           test_root_of<Root_of_2<CGAL::Quotient<CGAL::MP_Float> > >();

#ifdef CGAL_USE_GMP
  std::cout << "Testing Root_of_2<Gmpz>" << std::endl;
  result = result && test_root_of<Root_of_2<CGAL::Gmpz> >();

  std::cout << "Testing Root_of_2<Gmpq>" << std::endl;
  result = result && test_root_of<Root_of_2<CGAL::Gmpq> >();
#endif

#ifdef CGAL_USE_GMPXX
  std::cout << "Testing Root_of_2<mpq_class>" << std::endl;
  result = result && test_root_of<Root_of_2<mpq_class> >();

  std::cout << "Testing Root_of_2<Quotient<mpz_class> >" << std::endl;
  result = result && test_root_of<Root_of_2<CGAL::Quotient<mpz_class> > >();

  std::cout << "Testing Root_of_2<mpz_class>" << std::endl;
  result = result && test_root_of<Root_of_2<mpz_class> >();
#endif

#ifdef CGAL_USE_LEDA
  std::cout << "Testing Root_of_2<leda_real>" << std::endl;
  result = result && test_root_of<Root_of_2<leda_real> >();
#endif

  if (result) {
    std::cout << "OK" << std::endl;
    return 0;
  } else {
    std::cout << "ERROR" << std::endl;
    return -1;
  }
}
