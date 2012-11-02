/****************************************************************************
 * Core Library Version 1.7, August 2004
 * Copyright (c) 1995-2004 Exact Computation Project
 * All rights reserved.
 *
 * This file is part of CORE (http://cs.nyu.edu/exact/core/).
 * You can redistribute it and/or modify it under the terms of the GNU
 * General Public License as published by the Free Software Foundation,
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
 * File: extLong.h
 * Synopsis: 
 * 		An extended class for long
 *
 * Written by 
 *       Koji Ouchi <ouchi@simulation.nyu.edu>
 *       Chee Yap <yap@cs.nyu.edu>
 *       Igor Pechtchanski <pechtcha@cs.nyu.edu>,
 *       Vijay Karamcheti <vijayk@cs.nyu.edu>,
 *       Chen Li <chenli@cs.nyu.edu>
 *       Zilin Du <zilin@cs.nyu.edu>
 *
 * WWW URL: http://cs.nyu.edu/exact/
 * Email: exact@cs.nyu.edu
 *
 * $URL$
 * $Id$
 ***************************************************************************/

#ifndef _CORE_EXTLONG_H_
#define _CORE_EXTLONG_H_

#include <CGAL/CORE/Impl.h>
#include <CGAL/CORE/CoreAux.h>

namespace CORE { 

#ifndef LONG_MAX
#error "haven't define LONG_MAX"
#endif

#ifndef LONG_MIN
#error "haven't define LONG_MIN"
#endif

// LONG_MAX and LONG_MIN is assumed in this class:

const long EXTLONG_MAX = LONG_MAX;
const long EXTLONG_MIN = LONG_MIN + 1;
const long EXTLONG_NAN = LONG_MIN;
const unsigned long U_EXTLONG_MAX = LONG_MAX;

/// \class extLong
/// \brief extended long integer
class CGAL_CORE_EXPORT extLong {
private:
  long val;  ///< internal representation
  int  flag; ///< flags
  /**<  0 -- Normal;
        1 -- Overflow (positive);
       -1 -- Overflow (negative);
        2 -- NaN (sign can not be determined) */

  static void add(extLong& z, long x, long y);

public:

  /// \name Constructors
  //@{
  /// default constructor
  extLong();
  /// constructor for \c bool
  extLong(bool isNaN);
  /// constructor for \c int
  extLong(int);
  /// constructor for \c unsigned int
  extLong(unsigned int);
  /// constructor for \c long
  extLong(long);
  /// constructor for \c unsigned long
  extLong(unsigned long);
  //@}

  /// \name Arithmetic and assignment operators
  //@{
  extLong& operator +=(const extLong&);
  extLong& operator -=(const extLong&);
  extLong& operator *=(const extLong&);
  extLong& operator /=(const extLong&);
  //@}

  /// \name Incremental, Decremental, Unary minus operators
  //@{
  extLong& operator++();
  extLong  operator++(int);
  extLong& operator--();
  extLong  operator--(int);
  extLong  operator-() const;
  //@}

  /// \name Conversion Function
  //@{
  std::string toString() const {
    std::stringstream st;
    st << (*this);
    return st.str();
  }    
  long toLong() const;
  //@}

  /// \name Builtin functions
  //@{
  long asLong() const;
  bool isInfty() const;
  bool isTiny() const;
  bool isNaN() const;
  int  sign() const;
  /// comparison
  int cmp(const extLong &) const;
  //@}

  /// \name I/O Stream
  ///@{
  friend  CGAL_CORE_EXPORT std::ostream& operator <<(std::ostream&, const extLong&);
  //@}

  static const extLong& getNaNLong();
  static const extLong& getPosInfty();
  static const extLong& getNegInfty();
};



// constants (Globally)
#define CORE_NaNLong extLong::getNaNLong()
#define CORE_posInfty extLong::getPosInfty()
#define CORE_negInfty extLong::getNegInfty()

const extLong EXTLONG_ZERO(0);
const extLong EXTLONG_ONE(1);
const extLong EXTLONG_TWO(2);
const extLong EXTLONG_THREE(3);
const extLong EXTLONG_FOUR(4);
const extLong EXTLONG_FIVE(5);
const extLong EXTLONG_SIX(6);
const extLong EXTLONG_SEVEN(7);
const extLong EXTLONG_EIGHT(8);

// inline functions

//  private comparison function
inline int extLong::cmp(const extLong& x) const {
  if (isNaN() || x.isNaN()) {
    core_error("Two extLong NaN's cannot be compared!",
               __FILE__, __LINE__, false);
  }
  return (val == x.val) ? 0 : ((val > x.val) ? 1 : -1);
}

// default constructor (cheapest one)
inline extLong::extLong() : val(0), flag(0) {}

inline extLong::extLong(int i) : val(i), flag(0) {
  if (val == EXTLONG_MAX)
    flag = 1;
  else if (val <= EXTLONG_MIN)
    flag = -1;
}

inline extLong::extLong(unsigned int ui) : val(ui), flag(0) {
  if (val >= EXTLONG_MAX) {
    val  = EXTLONG_MAX;
    flag = 1;
  }
}

inline extLong::extLong(long l) : val(l), flag(0) {
  if (val >= EXTLONG_MAX)
    flag = 1;
  else if (val <= EXTLONG_MIN)
    flag = -1;
}

inline extLong::extLong(unsigned long u) {
  if (u >= U_EXTLONG_MAX) {
    val  = EXTLONG_MAX;
    flag = 1;
  } else {
    val = static_cast<long>(u);
    flag = 0;
  }
}

// isNaN defaults to false
inline extLong::extLong(bool isNaN) : val(0), flag(0) {
  if (isNaN) {
    val = EXTLONG_NAN;
    flag = 2;
  }
}

// comparison operators
inline bool operator== (const extLong& x, const extLong& y) {
  return x.cmp(y) == 0;
}

inline bool operator!= (const extLong& x, const extLong& y) {
  return x.cmp(y) != 0;
}

inline bool operator< (const extLong& x, const extLong& y) {
  return x.cmp(y) < 0;
}

inline bool operator<= (const extLong& x, const extLong& y) {
  return x.cmp(y) <= 0;
}

inline bool operator> (const extLong& x, const extLong& y) {
  return x.cmp(y) > 0;
}

inline bool operator>= (const extLong& x, const extLong& y) {
  return x.cmp(y) >= 0;
}

//  arithmetic operators
inline extLong operator+ (const extLong& x, const extLong& y) {
  return extLong(x)+=y;
}

inline extLong operator- (const extLong& x, const extLong& y) {
  return extLong(x)-=y;
}

inline extLong operator* (const extLong& x, const extLong& y) {
  return extLong(x)*=y;
}

inline extLong operator/ (const extLong& x, const extLong& y) {
  return extLong(x)/=y;
}

inline extLong& extLong::operator++ () {
  *this += 1;
  return *this;
}

inline extLong extLong::operator++ (int) {
  extLong r(*this);
  *this += 1;
  return r;
}

inline extLong& extLong::operator-- () {
  *this -= 1;
  return *this;
}

inline extLong extLong::operator-- (int) {
  extLong r(*this);
  *this -= 1;
  return r;
}

//  conversion to long
inline long extLong::toLong() const {
  return val;
}

// builtin functions
inline long extLong::asLong() const {
  return val;
}

inline bool extLong::isInfty() const {
  return (flag == 1);
}

inline bool extLong::isTiny() const {
  return (flag == -1);
}

inline bool extLong::isNaN() const {
  return (flag == 2);
}

} //namespace CORE
#endif // _CORE_EXTLONG_H_
