/******************************************************************
 * Core Library Version 1.5, August 2002
 * Copyright (c) 1995-2002 Exact Computation Project
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
 * $Id$
 *****************************************************************/

#include "extLong.h"

CORE_BEGIN_NAMESPACE

#ifndef CORE_ENABLE_INLINES
#include "extLong.inl"
#endif

const long nearHalfInfty    = LONG_MAX / 2;
const long farHalfInfty     = LONG_MAX - nearHalfInfty;
const long nearHalfTinyLong = (LONG_MIN + 1) / 2;
const long farHalfTinyLong  = LONG_MIN + 1 - nearHalfTinyLong;

const extLong& extLong::getNaNLong() {
  static extLong NaNLong(true);
  return NaNLong;
}

const extLong& extLong::getPosInfty() {
  static extLong posInfty(LONG_MAX);
  return posInfty;
}

const extLong& extLong::getNegInfty() {
  static extLong negInfty(LONG_MIN+1);
  return negInfty;
}

long extLong::add4Long(long x, long y) {
  long z;
  
  if (x >= farHalfInfty) {//  x is large
    long dx = x - farHalfInfty;
    
    if (y >= nearHalfInfty - dx)
      z = LONG_MAX; 
    else
      z = x + y;
  } else if (x <= farHalfTinyLong) {//  x is small
    long dx = farHalfTinyLong - x;
    
    if (y <= nearHalfTinyLong + dx)
      z = LONG_MIN + 1; 
    else
      z = x + y;
  } else { //  x is not small nor large
    if (y >= farHalfInfty) {//  y is large
      long dy = y - farHalfInfty;
      
      if (x >= nearHalfInfty - dy)
        z = LONG_MAX;
      else
        z = x + y;
    } else if (y <= farHalfTinyLong) {//  y is small
      long dy = farHalfTinyLong - y;
      
      if (x <= nearHalfTinyLong + dy)
        z = LONG_MIN + 1;
      else
        z = x + y;
    } else {//  both x and y are not small nor large
      z = x + y;
    }
  }
  return z;
}

long extLong::sub4Long(long x, long y) {
  long z;
  
  if (x >= farHalfInfty) {//  x is large
    long dx = x - farHalfInfty;
    
    if (y <= - nearHalfInfty + dx)
      z = LONG_MAX; 
    else
      z = x - y;
  } else if (x <= farHalfTinyLong) {//  x is small
    long dx = farHalfTinyLong - x;
    
    if (y >= - nearHalfTinyLong - dx)
      z = LONG_MIN + 1;
    else
      z = x - y;
  } else {//  x is not small nor large
    if (y <= - farHalfInfty) {//  y is small
      long dy = - farHalfInfty - y;
      
      if (x >= nearHalfInfty - dy)
        z = LONG_MAX; 
      else
        z = x - y;
    } else if (y >= - farHalfTinyLong) {//  y is large
      long dy = y + farHalfTinyLong;
      
      if (x <= nearHalfTinyLong + dy)
        z = LONG_MIN + 1; 
      else
        z = x - y;
    } else {//  both x and y are not small nor large
      z = x - y;
    }
  }
  return z;
}

//  arithmetic and assignment operators
extLong& extLong::operator+= (const extLong& y) {
  if (isNaN() || y.isNaN() || (flag * y.flag < 0)) {
#ifdef DEBUG
    if (flag * y.flag < 0) // want a message at the first creation of NaN
      core_error("extLong NaN Error in addition.", __FILE__, __LINE__, false);
#endif
    *this = CORE_NaNLong;
  } else if (flag == 1) {
    *this = CORE_posInfty;
  } else if (flag == -1) {
    *this  = CORE_negInfty;
  } else {
    // x has normal value now.
    if (y.flag == 1) { // y is positive infty
      *this = CORE_posInfty;
    } else if (y.flag == -1) { // y is negative infty
      *this = CORE_negInfty;
    } else { // y is normal
      val = add4Long(val, y.val);
      if (val == LONG_MAX)
        flag = 1;
      else if (val <= LONG_MIN + 1)
        flag = - 1;
      else
        flag = 0;
    }
  }
  return *this;
}

extLong& extLong::operator-= (const extLong& y) {
  if (isNaN() || y.isNaN() || (flag * y.flag > 0)) {
#ifdef DEBUG
    if (flag * y.flag > 0)	// want a message at the first creation of NaN
     core_error("extLong NaN Error in subtraction.", __FILE__, __LINE__, false);
#endif
    *this = CORE_NaNLong;
  } else if (flag == 1) {
    *this  = CORE_posInfty;
  } else if (flag == -1) {
    *this  = CORE_negInfty;
  } else {
    // x is normal now
    if (y.flag == 1) { // y is positive infty
      *this = CORE_negInfty;
    } else if (y.flag == -1) { // y is negative infty
      *this = CORE_posInfty;
    } else { // y is also normal.
      val = sub4Long(val, y.val);
      if (val == LONG_MAX)
        flag = 1;
      else if (val <= LONG_MIN + 1)
        flag = - 1;
      else
        flag = 0;
    }
  }
  return *this;
}

extLong& extLong::operator*= (const extLong& y) {
  if (isNaN() || y.isNaN()){
    *this = CORE_NaNLong;
  } else if ((flag != 0) || (y.flag != 0)) {
    int s = sign() * y.sign();
    switch (s) {
    case 0:
      val = 0; flag = 0;
      break;
    case 1: 
      *this = CORE_posInfty;
      break;
    case -1:
      *this = CORE_negInfty;
      break;
    case 2: 
    case -2: 		// Oct 30, 2002: this was a missing case
      *this = CORE_NaNLong;
      break;
    }
  } else { // flag == 0 and y.flag == 0
    double d = double(val) * double(y.val);
    long   p = val * y.val;
    if (fabs(d - p) <= fabs(d) * relEps) {
      val = p; flag = 0;
    } else if (d > LONG_MAX) {
      *this = CORE_posInfty;
    } else if (d < LONG_MIN + 1) {
      *this = CORE_negInfty;
    } else {
      core_error("extLong NaN Error in multiplication.", __FILE__, __LINE__, false);
      *this = CORE_NaNLong;
    }
  }
  return *this;
}

extLong& extLong::operator/= (const extLong& y) {
  if (isNaN() || y.isNaN() || ((flag != 0) && (y.flag != 0)) || (y.val == 0)) {
    if (y.val == 0)
     core_error("extLong NaN Error, Divide by Zero.", __FILE__, __LINE__, false);
    if ((flag !=0) && (y.flag !=0))
     core_error("extLong NaN Error, +/-Infinity/Infinity .", __FILE__,
		     	__LINE__, false);
    *this = CORE_NaNLong;
  } else if (flag != 0) { // y.flag == 0 now and y != 0
    int s = sign() * y.sign();
    switch (s) {
    case 0:
      val = 0; flag = 0;
      break;
    case 1: 
      *this = CORE_posInfty;
      break;
    case -1:
      *this = CORE_negInfty;
      break;
    case 2:  
    case -2:  			// Oct 30, 2002: this case was missing before
      *this = CORE_NaNLong;
      break;
    }
  } else { // flag == 0 and y.flag == 0
    val /= y.val; // no overflow in divisions
    flag = 0;
  }
  return *this;
}

//  unary minus
extLong extLong::operator- () const{
  extLong r;
  if (flag == 0) {
    r.val = - val; r.flag = 0;
  } else if (flag == 1) {
    r = CORE_negInfty;
  } else if (flag == -1) {
    r = CORE_posInfty;
  } else { // NaN
    r = CORE_NaNLong;
  }
  return r;
}

// sign
//    You should check "flag" before calling this, otherwise
//	you cannot interprete the returned value!
int extLong::sign() const {
  switch (flag) {
  case 0: 
    if (val > 0) return 1;
    else if (val < 0) return -1;
    else return 0;
    break;
  case 1:
    return 1;  // +infinity
    break;
  case -1:
    return -1;  // -infinity
    break;
  case 2:
    core_error("NaN Sign can not be determined!", __FILE__, __LINE__, false);
    return 2;
    break;
  default:
    core_error("Error sign flag!", __FILE__, __LINE__, false);
    return -2;
    break;
  }
}

//  stream operators
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

CORE_END_NAMESPACE
