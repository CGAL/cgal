// Copyright (c) 2002, 2003  Utrecht University (The Netherlands),
// ETH Zurich (Switzerland), Freie Universitaet Berlin (Germany),
// INRIA Sophia-Antipolis (France), Martin-Luther-University Halle-Wittenberg
// (Germany), Max-Planck-Institute Saarbruecken (Germany), RISC Linz (Austria),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Michael Hoffmann, Sylvain Pion

#ifndef CGAL_TEST_UTILITIES_H
#define CGAL_TEST_UTILITIES_H

#include <CGAL/basic.h>
#include <cassert>

CGAL_BEGIN_NAMESPACE

// Used to shut down some warnings about unused variables.
template < class T >
void use(const T&) {}

template < class NT >
bool
test_sqrt(const NT&, CGAL::Tag_true)
{
  NT sixteen(16);
  NT four(4);
  std::cout << "  sqrt()" << std::endl;
  if (CGAL_NTS sqrt(sixteen) != four) return false;
  CGAL::Sqrt<NT> s;
  if (s(sixteen) != four) return false;

  return true;
}

template < class NT >
bool
test_sqrt(const NT&, CGAL::Tag_false)
{
  return true;
}

template < class NT >
bool
test_gcd(const NT&, CGAL::Tag_true)
{
  // div
  NT eleven(11);
  NT three(3);
  std::cout << "  div()" << std::endl;
  if (CGAL_NTS div(eleven, three) != three) return false;
  CGAL::Div<NT> d;
  if (d(eleven, three) != three) return false;

  // gcd
  NT x(2737);
  NT y(5083);
  std::cout << "  gcd()" << std::endl;
  CGAL::Gcd<NT> gc;
  if (CGAL_NTS gcd(x, y) != NT(391)) return false;
  if (gc(x, y) != NT(391)) return false;

  typedef typename CGAL::Number_type_traits<NT>::Has_sqrt Has_sqrt;
  Has_sqrt has_sqrt = Has_sqrt();
  return test_sqrt(x, has_sqrt);
}

template < class NT >
bool
test_gcd(const NT &x, CGAL::Tag_false)
{
  typedef typename CGAL::Number_type_traits<NT>::Has_sqrt Has_sqrt;
  Has_sqrt has_sqrt = Has_sqrt();
  return test_sqrt(x, has_sqrt);
}


template < class NT >
bool
test_basic_operators(const NT&)
{
  std::cout << "  basic operators" << std::endl;
  NT zero(0);
  NT one(1);

  NT a = zero + one;
  NT b = zero - one;
  a = zero * one;
  a += b;
  a -= b;
  a *= b;
  a += 1;
  a -= 1;
  a *= 1;
  a = -b;

  bool d;
  d = a<b;
  d = a>b;
  d = a<=b;
  d = a>=b;
  d = a==b;
  d = a!=b;

  use(a); use(b); use(d);

  return true;
}

template < class NT >
bool
test_mixed_operators(const NT&)
{
  std::cout << "  mixed operators" << std::endl;

  NT zero = 0;
  NT one  = 1;

  // Comparison operators.
  if (   zero != 0  ) return false;
  if (   0 != zero  ) return false;
  if (! (zero == 0) ) return false;
  if (! (0 == zero) ) return false;

  if (   zero < 0   ) return false;
  if (   0 < zero   ) return false;
  if (! (zero >= 0) ) return false;
  if (! (0 >= zero) ) return false;

  if (! (zero <= 0) ) return false;
  if (! (0 <= zero) ) return false;
  if (   zero > 0   ) return false;
  if (   0 > zero   ) return false;

  // More tests for the semantics...
  if (! (zero != 1) ) return false;
  if (! (1 != zero) ) return false;
  if (   zero == 1  ) return false;
  if (   1 == zero  ) return false;

  if (! (zero < 1)  ) return false;
  if (   1 < zero   ) return false;
  if (   zero >= 1  ) return false;
  if (! (1 >= zero) ) return false;

  if (! (zero <= 1) ) return false;
  if (   1 <= zero  ) return false;
  if (   zero > 1   ) return false;
  if (! (1 > zero)  ) return false;

  // compare()
  if (CGAL_NTS compare(zero, 1)        != CGAL::SMALLER) return false;
  if (CGAL_NTS compare(1, zero)        != CGAL::LARGER)  return false;
  // Same, for expression templates
  if (CGAL_NTS compare(zero + zero, 1) != CGAL::SMALLER) return false;
  if (CGAL_NTS compare(1, zero + zero) != CGAL::LARGER)  return false;

  // The other operators +, -, *, /.
  if (zero + 1 != 1 ) return false;
  if (1 + zero != 1 ) return false;
  if (zero - 1 !=-1 ) return false;
  if (1 - zero != 1 ) return false;
  if (zero * 1 != 0 ) return false;
  if (1 * zero != 0 ) return false;
  if (zero / 1 != 0 ) return false;
  if (1 / one  != 1 ) return false;

  return true;
}

template < class NT >
bool
test_utilities(const NT& x)
{
  typedef typename CGAL::Number_type_traits<NT>::Has_division   Has_division;
  typedef typename CGAL::Number_type_traits<NT>::Has_sqrt       Has_sqrt;
  typedef typename CGAL::Number_type_traits<NT>::Has_gcd        Has_gcd;
  typedef typename CGAL::Number_type_traits<NT>::Has_exact_sqrt Has_exact_sqrt;
  typedef typename CGAL::Number_type_traits<NT>::Has_exact_division
                                                    Has_exact_division;
  typedef typename CGAL::Number_type_traits<NT>::Has_exact_ring_operations
                                                    Has_exact_ring_operations;


  NT zero(0);
  NT one(1);
  NT mone(-one);

  // is_zero
  std::cout << "  is_zero()" << std::endl;
  if (!CGAL_NTS is_zero(zero)) return false;
  CGAL::Is_zero<NT> iz;
  if (!iz(zero)) return false;

  // is_one
  std::cout << "  is_one()" << std::endl;
  if (!CGAL_NTS is_one(one)) return false;
  CGAL::Is_one<NT> io;
  if (!io(one)) return false;

  // is_positive
  std::cout << "  is_positive()" << std::endl;
  if (!CGAL_NTS is_positive(one)) return false;
  CGAL::Is_positive<NT> ip;
  if (!ip(one)) return false;

  // is_negative
  std::cout << "  is_negative()" << std::endl;
  if (CGAL_NTS is_negative(one)) return false;
  CGAL::Is_negative<NT> in;
  if (in(one)) return false;

  // sign
  std::cout << "  sign()" << std::endl;
  CGAL::Sgn<NT> sg;
  if (CGAL_NTS sign(one) != CGAL::POSITIVE) return false;
  if (sg(one) != CGAL::POSITIVE) return false;
  if (CGAL_NTS sign(zero) != CGAL::ZERO) return false;
  if (sg(zero) != CGAL::ZERO) return false;
  if (CGAL_NTS sign(mone) != CGAL::NEGATIVE &&
      CGAL_NTS sign(mone) != CGAL::POSITIVE) return false; // unsigned types
  if (sg(mone) != CGAL::NEGATIVE &&
      sg(mone) != CGAL::POSITIVE) return false; // unsigned types

  // abs
  std::cout << "  abs()" << std::endl;
  CGAL::Abs<NT> ab;
  if (CGAL_NTS abs(mone) != one &&
      CGAL_NTS abs(mone) != mone) // unsigned types :-)
    return false;
  if (ab(mone) != one && ab(mone) != mone)
    return false;
  if (ab(mone) != CGAL_NTS abs(mone) || ab(one) != CGAL_NTS abs(one))
    return false;

  // compare
  std::cout << "  compare()" << std::endl;
  CGAL::Compare<NT> co;
  if (CGAL_NTS compare(one, zero) != CGAL::LARGER) return false;
  if (CGAL_NTS compare(one, one) != CGAL::EQUAL) return false;
  if (CGAL_NTS compare(zero, one) != CGAL::SMALLER) return false;
  if (co(one, zero) != CGAL::LARGER) return false;
  if (co(one, one) != CGAL::EQUAL) return false;
  if (co(zero, one) != CGAL::SMALLER) return false;

  // to_double
  std::cout << "  to_double()" << std::endl;
  if (CGAL::to_double(zero) != 0.0) return false;
  if (CGAL::to_double(one)  != 1.0) return false;
  CGAL::To_double<NT> tdb;
  if (tdb(one) != 1.0) return false;

  // to_interval
  std::pair<double, double> I = CGAL::to_interval( NT(2) );
  if ( I.first  > 2.0) return false;
  if ( I.second < 2.0) return false;
  CGAL::To_interval<NT> tint;
  if (tint(one).first  > 1.0) return false;
  if (tint(one).second < 1.0) return false;

  // basic operators +,-,...
  if (!test_basic_operators(zero)) return false;

  // mixed operators +,-,...
  if (!test_mixed_operators(zero)) return false;

  // is_finite, is_valid
  std::cout << "  is_finite()" << std::endl;
  if (! CGAL::is_finite(zero)) return false;
  std::cout << "  is_valid()" << std::endl;
  if (! CGAL::is_valid(zero)) return false;

  // square
  std::cout << "  square()" << std::endl;
  NT two(2);
  NT four(4);
  CGAL::Square<NT> sq;
  if (CGAL_NTS square(two) != four) return false;
  if (sq(two) != four) return false;

  // Tests for number types that use the expression template mechanism
  // (e.g. mpz_class)
  if (!CGAL_NTS is_zero(zero+zero)) return false;
  if (!CGAL_NTS is_one(one-zero)) return false;
  if (!CGAL_NTS is_positive(one+zero)) return false;
  if (CGAL_NTS is_negative(one+one)) return false;
  if (CGAL_NTS sign(mone+mone) != CGAL::NEGATIVE) return false;
  if (CGAL_NTS abs(mone+mone) != one+one) return false;
  if (CGAL_NTS compare(zero+zero, one-zero) != CGAL::SMALLER) return false;
  if (CGAL_NTS compare(zero, one-zero) != CGAL::SMALLER) return false;
  if (CGAL_NTS compare(zero+zero, one) != CGAL::SMALLER) return false;
  if (CGAL::to_double(zero+zero) != 0.0) return false;
  if (! CGAL::is_finite(zero+zero)) return false;
  if (! CGAL::is_valid(zero+zero)) return false;
  if (CGAL_NTS square(two+two) != NT(16)) return false;

  // gcd
  typedef typename CGAL::Number_type_traits<NT>::Has_gcd Has_gcd;
  Has_gcd has_gcd = Has_gcd();
  return test_gcd(x, has_gcd);
}

CGAL_END_NAMESPACE

#endif
