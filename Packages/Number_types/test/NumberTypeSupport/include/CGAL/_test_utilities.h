// ============================================================================
//
// Copyright (c) 2002 The CGAL Consortium
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
// file          : utilites.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Michael Hoffmann
//
// coordinator   : MPI, Saarbruecken
// ============================================================================
 
#ifndef CGAL_TEST_UTILITIES_H
#define CGAL_TEST_UTILITIES_H

#include <CGAL/basic.h>
#include <cassert>

CGAL_BEGIN_NAMESPACE

template < class NT >
bool
test_sqrt(NT, CGAL::Tag_true)
{
  NT sixteen(16);
  NT four(4);
  std::cerr << "  sqrt()" << std::endl;
  if (CGAL_NTS sqrt(sixteen) != four) return false;
  CGAL::Sqrt<NT> s;
  if (s(sixteen) != four) return false;

  return true;
}

template < class NT >
bool
test_sqrt(NT, CGAL::Tag_false)
{
  return true;
}

template < class NT >
bool
test_gcd(NT z, CGAL::Tag_true)
{
  // div
  NT eleven(11);
  NT three(3);
  std::cerr << "  div()" << std::endl;
  if (CGAL_NTS div(eleven, three) != three) return false;
  CGAL::Div<NT> d;
  if (d(eleven, three) != three) return false;
  
  // gcd
  NT x(2737);
  NT y(5083);
  std::cerr << "  gcd()" << std::endl;
  CGAL::Gcd<NT> gc;
  if (CGAL_NTS gcd(x, y) != NT(391)) return false;
  if (gc(x, y) != NT(391)) return false;

  typedef typename CGAL::Number_type_traits<NT>::Has_sqrt Has_sqrt;
  Has_sqrt has_sqrt = Has_sqrt();
  return test_sqrt(x, has_sqrt);
}

template < class NT >
bool
test_gcd(NT x, CGAL::Tag_false)
{
  typedef typename CGAL::Number_type_traits<NT>::Has_sqrt Has_sqrt;
  Has_sqrt has_sqrt = Has_sqrt();
  return test_sqrt(x, has_sqrt);
}

template < class NT >
bool
test_utilities(NT x)
{
  NT zero(0);
  NT one(1);
  NT mone(-one);

  // is_zero
  std::cerr << "  is_zero()" << std::endl;
  if (!CGAL_NTS is_zero(zero)) return false;
  CGAL::Is_zero<NT> iz;
  if (!iz(zero)) return false;

  // is_one
  std::cerr << "  is_one()" << std::endl;
  if (!CGAL_NTS is_one(one)) return false;
  CGAL::Is_one<NT> io;
  if (!io(one)) return false;

  // is_positive
  std::cerr << "  is_positive()" << std::endl;
  if (!CGAL_NTS is_positive(one)) return false;
  CGAL::Is_positive<NT> ip;
  if (!ip(one)) return false;

  // is_negative
  std::cerr << "  is_negative()" << std::endl;
  if (CGAL_NTS is_negative(one)) return false;
  CGAL::Is_negative<NT> in;
  if (in(one)) return false;

  // sign
  std::cerr << "  sign()" << std::endl;
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
  std::cerr << "  abs()" << std::endl;
  CGAL::Abs<NT> ab;
  if (CGAL_NTS abs(mone) != one &&
      CGAL_NTS abs(mone) != mone) // unsigned types :-)
    return false;
  if (ab(mone) != one && ab(mone) != mone) 
    return false;
  if (ab(mone) != CGAL_NTS abs(mone) || ab(one) != CGAL_NTS abs(one))
    return false;

  // compare
  std::cerr << "  compare()" << std::endl;
  CGAL::Compare<NT> co;
  if (CGAL_NTS compare(one, zero) != CGAL::LARGER) return false;
  if (CGAL_NTS compare(one, one) != CGAL::EQUAL) return false;
  if (CGAL_NTS compare(zero, one) != CGAL::SMALLER) return false;
  if (co(one, zero) != CGAL::LARGER) return false;
  if (co(one, one) != CGAL::EQUAL) return false;
  if (co(zero, one) != CGAL::SMALLER) return false;

  // square
  std::cerr << "  square()" << std::endl;
  NT two(2);
  NT four(4);
  CGAL::Square<NT> sq;
  if (CGAL_NTS square(two) != four) return false;
  if (sq(two) != four) return false;

  typedef typename CGAL::Number_type_traits<NT>::Has_gcd Has_gcd;
  Has_gcd has_gcd = Has_gcd();
  return test_gcd(x, has_gcd);
}

CGAL_END_NAMESPACE

#endif
