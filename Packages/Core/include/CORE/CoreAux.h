/******************************************************************
 * Core Library Version 1.6, June 2003
 * Copyright (c) 1995-2002 Exact Computation Project
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
 * $Id$
 *****************************************************************/

#ifndef CORE_AUX_H
#define CORE_AUX_H

#include "CoreImpl.h"

CORE_BEGIN_NAMESPACE

#ifndef LONG_BIT // such as in Linux
  #define LONG_BIT (sizeof(long) * 8)
#endif

/// define machine epsilon for IEEE standard double
// const double eps = (ldexp(1.0, -53)); // fails to link on SunPro
#define CORE_EPS ((1.0/(1<<30))/(1<<23))
/// define relative machine epsilon for IEEE standard double
const double relEps = (1.0 + ldexp(1.0, -52));

/// template function returns the maximum value of two
template <class T>
CORE_INLINE const T& core_max(const T& a, const T& b)
{ return ((a > b) ? a : b); }

/// template function returns the minimum value of two
template <class T>
CORE_INLINE const T& core_min(const T& a, const T& b)
{ return ((a < b) ? a : b); }

/// template function returns the maximum value of three
template <class T>
CORE_INLINE const T& core_max(const T& a, const T& b, const T& c)
{ return ((a > b) ? core_max(a, c) : core_max(b, c)); }

/// template function swaps two values
template <class T>
CORE_INLINE  void core_swap(T& a, T& b)
{ T tmp; tmp = a;  a = b;  b = tmp; }

/// template function rotate three values
template <class T>
CORE_INLINE  void core_rotate(T& a, T& b, T& c)
{ T tmp; tmp = a;  a = b;  b = c; c = tmp; }

/// template function returns the minimum value of three
template <class T>
CORE_INLINE const T& core_min(const T& a, const T& b, const T& c)
{  return ((a < b) ? core_min(a, c) : core_min(b, c)); }

/// template function returns the absolute value
template <class T>
CORE_INLINE const T core_abs(const T& a)
{ return ((a < T(0)) ? -a : a); }

/// returns floor log base 2 of abs(x)
/**  CONVENTION: lg(0) = -1 */
CORE_INLINE int flrLg(long x);

/// returns floor log base 2 of unsigned long x
/**  CONVENTION: lg(0) = -1 */
CORE_INLINE int flrLg(unsigned long x);

/// returns ceiling log base 2 of abs(x)
/**  CONVENTION: lg(0) = -1 */
CORE_INLINE int clLg(long x);

/// returns ceiling log base 2 of unsigned long x
/**  CONVENTION: lg(0) = -1 */
CORE_INLINE int clLg(unsigned long x);

/// overloaded << to print out std::string
CORE_INLINE std::ostream& operator<< (std::ostream& o, const std::string& s);

/// print out an error message
CORE_INLINE void core_error(std::string msg,
		std::string filename, int lineno, bool ex);

/// gcd for machine type long
CORE_INLINE long gcd(long m, long n);

#ifdef CORE_ENABLE_INLINES
  #include "CoreAux.inl"
#endif

CORE_END_NAMESPACE
#endif
