/****************************************************************************
 * Core Library Version 1.7, August 2004
 * Copyright (c) 1995-2004 Exact Computation Project
 * All rights reserved.
 *
 * This file is part of CGAL (www.cgal.org).
 * You can redistribute it and/or modify it under the terms of the GNU
 * Lesser General Public License as published by the Free Software Foundation,
 * either version 3 of the License, or (at your option) any later version.
 *
 * Licensees holding a valid commercial license may use this file in
 * accordance with the commercial license agreement provided with the
 * software.
 *
 * This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
 * WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
 *
 *
 * File: extLong.cpp
 * Synopsis:
 *      The class extLong is basically a wrapper around the machine
 *      type long.  It is an important class to provide several
 *      additional facilities to detect overflows and undefined values.
 *      Future development includes extensions to level arithmetic
 *      (i.e., if a number overflows level i, we will go to level i+1).
 *      Level i representation of a number n is just i iterations
 *      of log_2 applied to n.
 *
 * Written by 
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
 ***************************************************************************/

#ifdef CGAL_HEADER_ONLY
#define CGAL_INLINE_FUNCTION inline
#else
#define CGAL_INLINE_FUNCTION
#endif

#include <CGAL/CORE/extLong.h>

namespace CORE { 

CGAL_INLINE_FUNCTION
const extLong& extLong::getNaNLong() {
  static const extLong NaNLong(true);
  return NaNLong;
}

CGAL_INLINE_FUNCTION
const extLong& extLong::getPosInfty() {
  static const extLong posInfty(EXTLONG_MAX);
  return posInfty;
}

CGAL_INLINE_FUNCTION
const extLong& extLong::getNegInfty() {
  static const extLong negInfty(EXTLONG_MIN);
  return negInfty;
}

CGAL_INLINE_FUNCTION
void extLong::add(extLong& z, long x, long y) {
  if (x > 0 && y > 0 && x >= EXTLONG_MAX - y) {
    z.val = EXTLONG_MAX;
    z.flag = 1;
  } else if (x < 0 && y < 0 && x <= EXTLONG_MIN - y) {
    z.val = EXTLONG_MIN;
    z.flag = -1;
  } else {
    z.val = x + y;
    z.flag = 0;
  }
}

//  arithmetic and assignment operators
CGAL_INLINE_FUNCTION
extLong& extLong::operator+= (const extLong& y) {
  if (flag == 2 || y.flag == 2 || (flag * y.flag < 0)) {
#ifdef CORE_DEBUG
    if (flag * y.flag < 0) //want a message at the first creation of NaN
      core_error("extLong NaN Error in addition.", __FILE__, __LINE__, false);
#endif

    *this = CORE_NaNLong;
  } else if (flag == 1 || y.flag == 1) { // one of them is +Inf
    *this = CORE_posInfty;
  } else if (flag == -1 || y.flag == -1) { // one of them is -Inf
    *this = CORE_negInfty;
  } else { // x and y are normal now
    add(*this, val, y.val);
  }
  return *this;
}

CGAL_INLINE_FUNCTION
extLong& extLong::operator-= (const extLong& y) {
  if (flag == 2 || y.flag == 2 || (flag * y.flag > 0)) {
#ifdef CORE_DEBUG
    if (flag * y.flag > 0) //want a message at the first creation of NaN
      core_error("extLong NaN Error in subtraction.", __FILE__, __LINE__, false);
#endif

    *this = CORE_NaNLong;
  } else if (flag == 1 || y.flag == -1) {
    *this = CORE_posInfty;
  } else if (flag == -1 || y.flag == 1) {
    *this = CORE_negInfty;
  } else {
    add(*this, val, -y.val);
  }
  return *this;
}

CGAL_INLINE_FUNCTION
extLong& extLong::operator*= (const extLong& y) {
  if (flag == 2 || y.flag == 2) {
    *this = CORE_NaNLong;
  } else if ((flag != 0) || (y.flag != 0)) {
    if (sign() * y.sign() > 0)
      *this = CORE_posInfty;
    else
      *this = CORE_negInfty;
  } else { // flag == 0 and y.flag == 0
    double d = double(val) * double(y.val);
    long   p = val * y.val;
    if (std::fabs(d - p) <= std::fabs(d) * relEps) {
      val = p;
      flag = 0;
    } else if (d > EXTLONG_MAX) {
      *this = CORE_posInfty;
    } else if (d < EXTLONG_MIN) {
      *this = CORE_negInfty;
    } else {
#ifdef CORE_DEBUG
      core_error("extLong NaN Error in multiplication.",__FILE__,__LINE__,false);
#endif
      *this = CORE_NaNLong;
    }
  }
  return *this;
}

CGAL_INLINE_FUNCTION
extLong& extLong::operator/= (const extLong& y) {
  if (flag==2 || y.flag==2 || ((flag != 0) && (y.flag != 0)) || (y.val == 0)) {
#ifdef CORE_DEBUG
    if (y.val == 0)
      core_error("extLong NaN Error, Divide by Zero.", __FILE__, __LINE__, false);
    else if ((flag !=0) && (y.flag !=0))
      core_error("extLong NaN Error, +/-Inf/Inf.", __FILE__, __LINE__, false);
#endif

    *this = CORE_NaNLong;
  } else if ((flag != 0) || (y.flag != 0)) { // y.flag == 0 now and y != 0
    if (sign() * y.sign() > 0)
      *this = CORE_posInfty;
    else
      *this = CORE_negInfty;
  } else { // flag == 0 and y.flag == 0
    val /= y.val; // no overflow in divisions
    flag = 0;
  }
  return *this;
}

//  unary minus
CGAL_INLINE_FUNCTION
extLong extLong::operator- () const {
  if (flag == 0)
    return extLong(-val);
  else if (flag == 1)
    return CORE_negInfty;
  else if (flag == -1)
    return CORE_posInfty;
  else // NaN
    return CORE_NaNLong;
}

// sign
//    You should check "flag" before calling this, otherwise
//    you cannot interprete the returned value!
CGAL_INLINE_FUNCTION
int extLong::sign() const {
  if (flag == 2)
    core_error("NaN Sign can not be determined!", __FILE__, __LINE__, false);
  return ((val == 0) ? 0 : ((val > 0) ? 1 : -1));
}

//  stream operators
CGAL_INLINE_FUNCTION
std::ostream& operator<< (std::ostream& o, const extLong& x) {
  if (x.flag == 1)
    o << " infty ";
  else if (x.flag == - 1)
    o << " tiny ";
  else if (x.flag == 2)
    o << " NaN ";
  else
    o << x.val;
  return o;
}

} //namespace CORE
