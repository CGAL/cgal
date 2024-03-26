/****************************************************************************
 * Core Library Version 1.7, August 2004
 * Copyright (c) 1995-2004 Exact Computation Project
 * All rights reserved.
 *
 * This file is part of CGAL (www.cgal.org).
 *
 * File: BigInt.h
 * Synopsis:
 *                 a wrapper class for mpz from GMP
 *
 * Written by
 *       Chee Yap <yap@cs.nyu.edu>
 *       Chen Li <chenli@cs.nyu.edu>
 *       Zilin Du <zilin@cs.nyu.edu>
 *
 * WWW URL: https://cs.nyu.edu/exact/
 * Email: exact@cs.nyu.edu
 *
 * $URL$
 * $Id$
 * SPDX-License-Identifier: LGPL-3.0-or-later
 ***************************************************************************/
#ifndef _CORE_BIGINT_H_
#define _CORE_BIGINT_H_

#include <CGAL/boost_mp_type.h>
#include <CGAL/CORE/RefCount.h>
#include <CGAL/CORE/MemoryPool.h>
#include <string>


#if !(defined(CGAL_CORE_USE_BOOST_BACKEND) && BOOST_VERSION > 107900 && defined(CGAL_USE_BOOST_MP)) && !defined(CGAL_DISABLE_GMP)
#define CGAL_CORE_USE_GMP_BACKEND 1
#endif

namespace CORE {

#ifdef CGAL_CORE_USE_GMP_BACKEND
  typedef boost::multiprecision::mpz_int BigInt;
#else
  typedef  boost::multiprecision::cpp_int BigInt;
#endif


inline int cmp(const BigInt& x, const BigInt& y) {
  return x.compare(y);
}


inline int set_str(BigInt& a, const char* s) {
    a = BigInt(s);
    return 0;  // should be -1 if not correct in the base (we ignore)
  }

  /// longValue
inline long longValue(const BigInt& a) {
  return a.convert_to<long>();
}

  /// doubleValue
inline double doubleValue(const BigInt& a) {
  return a.convert_to<double>();
}

  /// isEven
inline bool isEven(const BigInt& z) {
  return bit_test(z,0) == 0;
}
/// isOdd
inline bool isOdd(const BigInt& z) {
  return bit_test(z,0) == 1;
}

inline bool isDivisible(const BigInt& x, const BigInt& y) {
    BigInt q, r;
    divide_qr(x, y, q, r);
    return r.is_zero();
}

inline bool isDivisible(int x, int y) {
  return x % y == 0;
}

inline bool isDivisible(long x, long y) {
  return x % y == 0;

}
  /// get exponent of power 2
inline std::size_t getBinExpo(const BigInt& z) {
    if (z.is_zero()) {
        return (std::numeric_limits<std::size_t>::max)();
    }
    return lsb(abs(z));
}

  // bit length
inline std::size_t bitLength(const BigInt& a){
    if (a.is_zero()) {
        return 0;
    }
  return msb(abs(a))+1;
}

/// floorLg -- floor of log_2(a)
/** Convention: a=0, floorLg(a) returns -1.
 *  This makes sense for integer a.
 */
inline long floorLg(const BigInt& a) {
  assert(std::size_t((std::numeric_limits<long>::max)()) > bitLength(a));
  return (sign(a) == 0) ? (-1) : static_cast<long>(bitLength(a)-1);
}


/// div_rem
inline void div_rem(BigInt& q, BigInt& r, const BigInt& a, const BigInt& b) {
  divide_qr(a, b, q, r);
}


  /// ulongValue
inline unsigned long ulongValue(const BigInt& a) {
    assert(a >= BigInt(0));
  return a.convert_to<unsigned long>();
}

  /// exact div
inline void divexact(BigInt& z, const BigInt& x, const BigInt& y) {
  BigInt r;
  divide_qr(x, y, z, r );  // was void mpz_divexact (mpz_t q, const mpz_t n, const mpz_t d)   Is this faster?
  assert(r.is_zero());
}

// Chee (1/12/2004)   The definition of div_exact(x,y) next
//   ensure that in Polynomials<NT> works with both NT=BigInt and NT=int:
inline BigInt div_exact(const BigInt& x, const BigInt& y) {
  BigInt z;             // precodition: isDivisible(x,y)
  divexact(z, x, y); // z is set to x/y;
  return z;
}

inline int div_exact(int x, int y) {
  return x/y;  // precondition: isDivisible(x,y)
}

inline long div_exact(long x, long y) {
  return x/y;  // precondition: isDivisible(x,y)
}

inline BigInt gcd(const BigInt& a, const BigInt& b){
  return boost::multiprecision::gcd(a,b);
}

/// ceilLg -- ceiling of log_2(a) where a=BigInt, int or long
/** Convention: a=0, ceilLg(a) returns -1.
 *  This makes sense for integer a.
 */
inline long ceilLg(const BigInt& a) {
  if (sign(a) == 0)
    return -1;
  assert(std::size_t((std::numeric_limits<long>::max)()) > bitLength(a));
  std::size_t len = static_cast<long>(bitLength(a));

  return (lsb(abs(a)) == len - 1) ? (static_cast<long>(len) - 1) : static_cast<long>(len);
}

inline long ceilLg(long a) { // need this for Polynomial<long>
  return ceilLg(BigInt(a));
}

inline long ceilLg(int a) { // need this for Polynomial<int>
  return ceilLg(BigInt(a));
}


/// negate
inline void negate(BigInt& a) {
  a= - a;
}

/// get exponent of power k
inline void getKaryExpo(const BigInt& z, BigInt& m, int& e, unsigned long uk) {
    BigInt k(uk), q, r;
    e = 0;
    m = z;
    for(;;) {
        divide_qr(m, k, q, r);
        if (!r.is_zero()) break;
        m = q;
        ++e;
    }
}

inline void power(BigInt& c, const BigInt& a, unsigned long ul) {
  c = pow(a, ul);
}

} // namespace CORE



#endif // _CORE_BIGINT_H_
