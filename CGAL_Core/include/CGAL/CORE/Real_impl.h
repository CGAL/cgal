/****************************************************************************
 * Core Library Version 1.7, August 2004
 * Copyright (c) 1995-2004 Exact Computation Project
 * All rights reserved.
 *
 * This file is part of CGAL (www.cgal.org).
 *
 * File: Real.cpp
 * Synopsis: The Real class is a superclass for all the number
 *           systems in the Core Library (int, long, float, double,
 *           BigInt, BigRat, BigFloat, etc)
 *
 * Written by
 *       Koji Ouchi <ouchi@simulation.nyu.edu>
 *       Chee Yap <yap@cs.nyu.edu>
 *       Chen Li <chenli@cs.nyu.edu>
 *       Zilin Du <zilin@cs.nyu.edu>
 *       Sylvain Pion <pion@cs.nyu.edu>
 *
 * WWW URL: http://cs.nyu.edu/exact/
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

#include <CGAL/disable_warnings.h>

#include <ctype.h>
#include <CGAL/CORE/Real.h>
#include <CGAL/tss.h>
#ifdef CGAL_HEADER_ONLY
#include <CGAL/CORE/BigFloat.h> // for FiveTo
#endif

namespace CORE {

CGAL_INLINE_FUNCTION
const Real& Real::getZero() {
  init_CORE();
  CGAL_STATIC_THREAD_LOCAL_VARIABLE(Real, Zero, 0);
  return Zero;
}

CGAL_INLINE_FUNCTION
BigInt floor(const Real& r, Real &sub) {
  BigInt f = r.approx(CORE_INFTY, 2).BigIntValue();
  sub = r-f;
  // Adjustment
  if (sub<0)
    ++sub, --f;
  if (sub>=1)
    --sub, ++f;
  CGAL_assertion(sub >=0 && sub<1);
  return f;
}

// pow(r,n) and power(r, n) are the same function
//
CGAL_INLINE_FUNCTION
Real pow(const Real& r, unsigned long n) {
  if (n == 0)
    return Real(1);
  else if (n == 1)
    return r;
  else {
    Real x = r;
    while ((n % 2) == 0) { // n is even
      x *= x;
      n >>= 1;
    }
    Real u = x;
    while (true) {
      n >>= 1;
      if (n == 0)
        return u;
      x *= x;
      if ((n % 2) == 1) // n is odd
        u *= x;
    }
    //return u; // unreachable
  }
}//pow

extern BigInt FiveTo(unsigned long exp);

// Construct Real from String
// Note:
//         -- Zilin Du: 06/03/2003
//         -- Original it is the code for Real's constructor for "const char*".
//            I change it to a function so that two constrcutors can share the code.
//            now it is private and no default value.
//
//   --Default value of the argument "prec" is get_static_defInputDigits()
//   --If prec is CORE_posInfty, then the input is
//        read in exactly.  Otherwise, we convert to a RealBigFloat
//        with absolute error at most 10^{-prec}

// Constructor Real( char *str, extLong& prec)
//        is very similar to
//                BigFloatRep::fromString( char *str, extLong& prec);
// Differences:
//        In BigFloat(str, prec), the value of prec cannot be infinity, and
//                it defaults to defBigFloatInputDigits;
//        In Real(str, prec), the value of prec is allowed to be infinity, and
//                it defaults to defInputDigits.
//
// Why do we have the two versions?  It allows us to use the BigFloat class
//        directly, without relying on Real class.

CGAL_INLINE_FUNCTION
void Real::constructFromString(const char *str, const extLong& prec )
// NOTE: prec defaults to get_static_defInputDigits() (see Real.h)
{
  //        8/8/01, Chee and Zilin: add a new rational string format:
  //                this format is indicated by the presence of a slash "/"
  //                Moreover, the value of prec is ignored (basically
  //                assumed to be infinity).

  if (std::strchr(str, '/') != nullptr) {        // this is a rational number
    rep = new RealBigRat(BigRat(str));
    return;
  }

  const char *e = std::strchr(str, 'e');
  int dot = 0;
  long e10 = 0;
  if (e != nullptr)
    e10 = std::atol(e+1);        // e10 is decimal precision of the input string
  // i.e., input is A/10^{e10}.
  else {
    e = str + std::strlen(str);
#ifdef CORE_DEBUG
    CGAL_assertion(*e == '\0');
#endif
  }

  const char *p = str;
  if (*p == '-' || *p == '+')
    p++;
  BigInt m(0);

  for (; p < e; p++) {
    if (*p == '.') {
      dot = 1;
      continue;
    }
    m = m * 10 + (*p - '0');
    if (dot)
      e10--;
  }

  long t = (e10 < 0) ? -e10 : e10;
  BigInt one(1);
  BigInt ten = FiveTo(t) * (one << static_cast<unsigned long>(t));
  if (*str == '-')
    m = -m;
  if (e10 >= 0) {
    // convert exactly from integer numbers
    m *= ten;
    rep = new RealBigInt(m);
  } else { // e10 < 0,  fractional numbers
    // HERE IS WHERE WE USE THE SYSTEM CONSTANT
    //               get_static_defInputDigits()
    // Note: get_static_defInputDigits() should be at least log_2(10).
    //       We default get_static_defInputDigits() to 4.
    //std::cout << "(m,ten)=" << m << "," << ten << std::endl;
    BigRat r(m, ten);
    if (prec.isInfty()) { // convert exactly! to a big rational
      rep = new RealBigRat(r);
    } else {
      // convert approximately, to a BigFloat within the
      // specified precision:
      // BigFloat bf(r, CORE_posInfty, prec * lgTenM) ;
      BigFloat bf(r, CORE_posInfty, prec * 4) ;
      rep = new RealBigFloat(bf);
    }
  }
}// Real(str, prec)

// The operator >>(i,x) calls the constructor Real(char*)
CGAL_INLINE_FUNCTION
std::istream& operator >>(std::istream& i, Real& x) {
  int size = 20;
  char *str = new char[size];
  char *p = str;
  char c;
  int d = 0, e = 0, s = 0;
  //  int done = 0;

  // Chen Li: fixed a bug, the original statement is
  //  for (i.get(c); c == ' '; i.get(c));
  // use isspace instead of testing c == ' ', since it must also
  // skip tab, catridge/return, etc.
  // Change to:
  //  int status;
  do {
    i.get(c);
  } while (!i.eof() && isspace(c)); /* loop if met end-of-file, or
                             char read in is white-space. */
  // Chen Li,
  // original "if (c == EOF) ..." is unsafe since c is of char type and
  // EOF is of int tyep with a negative value -1

  if (i.eof()) {
    i.clear(std::ios::eofbit | std::ios::failbit);
    delete [] str;
    return i;
  }

  // the current content in "c" should be the first non-whitespace char
  if (c == '-' || c == '+') {
    *p++ = c;
    i.get(c);
  }

  for (; isdigit(c) || (!d && c=='.') ||
       (!e && c=='e') || (!s && (c=='-' || c=='+')); i.get(c)) {
    if (!i) break;
    if (!e && (c == '-' || c == '+'))
      break;
    // Chen Li: put one more rule to prohibite input like
    //  xxxx.xxxe+xxx.xxx:
    if (e && (c == '.'))
      break;
    if (p - str == size) {
      char *t = str;
      str = new char[size*2];
      std::memcpy(str, t, size);
      delete [] t;
      p = str + size;
      size *= 2;
    }
#ifdef CORE_DEBUG
    CGAL_assertion((p-str) < size);
#endif

    *p++ = c;
    if (c == '.')
      d = 1;
    // Chen Li: fix a bug -- the sign of exponent can not happen before
    // the character "e" appears! It must follow the "e' actually.
    //    if (e || c == '-' || c == '+') s = 1;
    if (e)
      s = 1;
    if (c == 'e')
      e = 1;
  }

  if (!i && !i.eof()) {
    delete [] str;
    return i;
  }
  // chenli: make sure that the p is still in the range
  if (p - str >= size) {
    std::ptrdiff_t len = p - str;
    char *t = str;
    str = new char[len + 1];
    std::memcpy(str, t, len);
    delete [] t;
    p = str + len;
  }

#ifdef CORE_DEBUG
  CGAL_assertion(p - str < size);
#endif

  *p = '\0';
  i.putback(c);
  i.clear();
  // old: x = Real(str, i.precision()); // use precision of input stream.
  x = Real(str);  // default precision = get_static_defInputDigits()
  delete [] str;
  return i;
}//operator >> (std::istream&, Real&)


} //namespace CORE

#include <CGAL/enable_warnings.h>
