/****************************************************************************
 * Core Library Version 1.7, August 2004
 * Copyright (c) 1995-2004 Exact Computation Project
 * All rights reserved.
 *
 * This file is part of CGAL (www.cgal.org).
 *
 * File: CoreAux.cpp
 * Synopsis:
 *       Auxiliary routines such as ceiling of log_2, etc.
 *       they are not specific to any Core classes.
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

#ifdef CGAL_HEADER_ONLY
#define CGAL_INLINE_FUNCTION inline
#else
#define CGAL_INLINE_FUNCTION
#endif

#include <CGAL/use.h>
#include <CGAL/CORE/CoreAux.h>

namespace CORE {

////////////////////////////////////////////////////////////
//  More useful functions to implement:
//
//  To convert digits into bits:
//      given X, compute X * log_2(10)
//  To convert bits into digits:
//      given X, compute X * log_10(2)
//
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// flrLg(x)
//      returns floor log base 2 of abs(x)
// CONVENTION lg(0) = -1        (Slight improvement, Zilin/Chee 8/5/02)
////////////////////////////////////////////////////////////
CGAL_INLINE_FUNCTION
int flrLg(long x) {
  if (x == LONG_MIN) {
    // special treatment as -LONG_MIN would be not representable as "long"
    return LONG_BIT - 1;
  } else {
    //  1 <= |x| <= LONG_MAX
    if (x < 0)
      x = - x;

    int lg = -1;
    while (x > 0) {
      lg++;
      x >>= 1;
    }
    return lg;
  }
}

////////////////////////////////////////////////////////////
// floor log base 2 of unsigned long x
// CONVENTION lg(0) = -1        (Slight improvement, Zilin/Chee 8/5/02)
////////////////////////////////////////////////////////////
CGAL_INLINE_FUNCTION
int flrLg(unsigned long x) {
  int lg = -1;
  while (x > 0) {
    lg++;
    x >>= 1;
  }
  return lg;
}

////////////////////////////////////////////////////////////
// ceiling log base 2 of abs(x)
// CONVENTION lg(0) = -1        (Slight improvement, Zilin/Chee 8/5/02)
////////////////////////////////////////////////////////////
CGAL_INLINE_FUNCTION
int clLg(long x) {
  if (x == LONG_MIN)
    return LONG_BIT - 1;
  if (x < 0)
    x = -x;                 // use absolute value
  if (x > (LONG_MAX >> 1))         // the leading data bit is 1
    return (LONG_BIT - 1);        // exclude the sign bit
  if (x >= 2)
    return flrLg((unsigned long)((x << 1) - 1));
  // SINCE ceilLog_2(x) = floorLog_2(2x-1) for x>=2
  if (x == 1)
    return 0;
  return -1;                        // x must be 0 here
}

////////////////////////////////////////////////////////////
// ceiling log base 2 of unsigned long x
// CONVENTION lg(0) = -1
////////////////////////////////////////////////////////////
CGAL_INLINE_FUNCTION
int clLg(unsigned long x) {
  if (x > (ULONG_MAX >> 1))        // the leading bit is 1
    return LONG_BIT;
  if (x >= 2)
    return flrLg((x << 1) - 1);
  // SINCE ceilLog_2(x) = floorLog_2(2x-1) for x>=2.
  if (x == 1)
    return 0;
  return -1;        // x must be equal to 0
}

/// gcd for machine type long
/** This is needed when we construct Polynomials with int or long coefficients */
CGAL_INLINE_FUNCTION
long gcd(long m, long n) {
  if (m == 0)
    return core_abs(n);
  if (n == 0)
    return core_abs(m);
  long p = core_abs(m);
  long q = core_abs(n);
  if (p<q)
    core_swap(p, q);
  while (q>0) {
    long r = p % q;
    p = q;
    q = r;
  }
  return p;
}


/// implements the "integer mantissa" function
//      (See CORE_PATH/progs/ieee/frexp.cpp for details)
CGAL_INLINE_FUNCTION
double IntMantissa(double d) {
        int e;
        return std::ldexp(std::frexp(d, &e), 53);
}

/// implements the "integer exponent" function
//      (See CORE_PATH/progs/ieee/frexp.cpp for details)
CGAL_INLINE_FUNCTION
int IntExponent(double d) {
        int e;
        std::frexp(d, &e);
        return e-53;
}


} //namespace CORE
