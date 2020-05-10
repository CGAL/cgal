/****************************************************************************
 * Core Library Version 1.7, August 2004
 * Copyright (c) 1995-2004 Exact Computation Project
 * All rights reserved.
 *
 * This file is part of CGAL (www.cgal.org).
 *
 * File: CoreAux.h
 * Synopsis:
 *      Auxilliary functions
 *
 * Written by
 *       Chee Yap <yap@cs.nyu.edu>
 *       Chen Li <chenli@cs.nyu.edu>
 *       Zilin Du <zilin@cs.nyu.edu>
 *
 * WWW URL: http://cs.nyu.edu/exact/
 * Email: exact@cs.nyu.edu
 *
 * $URL$
 * $Id$
 * SPDX-License-Identifier: LGPL-3.0-or-later
 ***************************************************************************/

#ifndef _CORE_COREAUX_H_
#define _CORE_COREAUX_H_

#include <iostream>
#include <fstream>
#include "CGAL/CORE/Impl.h"

namespace CORE {

#ifndef LONG_BIT // such as in Linux
  #define LONG_BIT (sizeof(long) * 8)
#endif

/// CORE_EPS is unit roundoff for IEEE standard double, i.e., 2^{-53}.
// NOTES:
// (1) CORE_EPS is used in our Floating Point Filter (Filter.h)
// (2) 2^{-53} is called "unit roundoff" and
//         is the roundoff error for numbers in the range [1,2).
//          "Machine epsilon" is 2*CORE_EPS = 2^{-52}.  It is the
//          smallest gap between two normal machine numbers  --Chee 8/2003
//
// const double eps = (ldexp(1.0, -53)); // fails to link on SunPro
#define CORE_EPS ((1.0/(1<<30))/(1<<23))
//
#define CORE_MACHINE_EPS ((1.0/(1<<30))/(1<<22))

/// relEps is relative error for IEEE standard double, 1+2^{-52}.
const double relEps = (1.0 + std::ldexp(1.0, -52));

/// CORE_DIAGFILE is used for all warning and error messages
const char* const CORE_DIAGFILE = "Core_Diagnostics";  // global file name

/// template function returns the maximum value of two
template <class T>
inline const T& core_max(const T& a, const T& b) {
  return ((a > b) ? a : b);
}

/// template function returns the minimum value of two
template <class T>
inline const T& core_min(const T& a, const T& b) {
  return ((a < b) ? a : b);
}

/// template function returns the maximum value of three
template <class T>
inline const T& core_max(const T& a, const T& b, const T& c) {
  return ((a > b) ? core_max(a, c) : core_max(b, c));
}

/// template function swaps two values
template <class T>
inline  void core_swap(T& a, T& b) {
  T tmp;
  tmp = a;
  a = b;
  b = tmp;
}

/// template function rotate three values
template <class T>
inline  void core_rotate(T& a, T& b, T& c) {
  T tmp;
  tmp = a;
  a = b;
  b = c;
  c = tmp;
}

/// template function returns the minimum value of three
template <class T>
inline const T& core_min(const T& a, const T& b, const T& c) {
  return ((a < b) ? core_min(a, c) : core_min(b, c));
}

/// template function returns the absolute value
template <class T>
inline const T core_abs(const T& a) {
  return ((a < T(0)) ? -a : a);
}

/// returns floor log base 2 of abs(x)
/**  CONVENTION: lg(0) = -1 */
CGAL_CORE_EXPORT int flrLg(long x);

/// returns floor log base 2 of unsigned long x
/**  CONVENTION: lg(0) = -1 */
CGAL_CORE_EXPORT int flrLg(unsigned long x);

/// returns ceiling log base 2 of abs(x)
/**  CONVENTION: lg(0) = -1 */
CGAL_CORE_EXPORT int clLg(long x);

/// returns ceiling log base 2 of unsigned long x
/**  CONVENTION: lg(0) = -1 */
CGAL_CORE_EXPORT int clLg(unsigned long x);

/// gcd for machine type long
CGAL_CORE_EXPORT long gcd(long m, long n);

/// abs for int type
inline int abs(int x) {
  return (x>=0) ? x : (-x);
}

/// abs for long type
inline long abs(long x) {
  return (x>=0) ? x : (-x);
}

/// sign for int type
inline int sign(int x) {
  return (x==0) ? 0 : ((x>0) ? 1 : (-1));
}

/// sign for long type
inline long sign(long x) {
  return (x==0) ? 0 : ((x>0) ? 1 : (-1));
}

/// overloaded << to print out std::string
inline std::ostream& operator<< (std::ostream& o, const std::string& s) {
  o << s.c_str();
  return o;
}

/// implements the "integer mantissa" function
//      (See CORE_PATH/progs/ieee/frexp.cpp for details)
CGAL_CORE_EXPORT double IntMantissa(double d);

/// implements the "integer exponent" function
//      (See CORE_PATH/progs/ieee/frexp.cpp for details)
CGAL_CORE_EXPORT int IntExponent(double d);

/// Writes out an error or warning message in the local file CORE_DIAGFILE
/** If last argument (err) is TRUE, then this is considered an error
 *  (not just warning).  In this case, the message is also printed in
 *  std::cerr, using std::perror().
 *  */
CGAL_CORE_EXPORT void core_error(std::string msg, std::string file, int lineno, bool err);

/// This is for debugging messages
inline void core_debug(std::string msg){
  std::cout << __FILE__ << "::" << __LINE__ << ": " << msg
            << std::endl;
}


} //namespace CORE

#ifdef CGAL_HEADER_ONLY
#include <CGAL/CORE/CoreAux_impl.h>
#endif // CGAL_HEADER_ONLY

#endif // _CORE_COREAUX_H_
