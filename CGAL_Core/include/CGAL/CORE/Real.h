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
 * File: Real.h
 * 
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
 ***************************************************************************/
#ifndef _CORE_REAL_H_
#define _CORE_REAL_H_
#include "RealRep.h"

namespace CORE { 
// class Real
typedef RCImpl<RealRep> RCReal;
class Real : public RCReal {
public:
  Real(int i=0) : RCReal(new RealLong(i)) {}
  Real(unsigned int ui) : RCReal(NULL) {
    (ui<=INT_MAX) ? (rep=new RealLong(static_cast<int>(ui))) : (rep=new RealBigInt(ui));
  }
  Real(long l) : RCReal(new RealLong(l)) {}
  Real(unsigned long ul) : RCReal(NULL) {
    (ul<=LONG_MAX) ? (rep=new RealLong(static_cast<long>(ul))) : (rep=new RealBigInt(ul));
  }
  Real(float f) : RCReal(new RealDouble(f)) {}
  Real(double d) : RCReal(new RealDouble(d)) {}
  Real(const BigInt& I) : RCReal(new RealBigInt(I)) {}
  Real(const BigRat& R) : RCReal(new RealBigRat(R)) {}
  Real(const BigFloat& F) : RCReal(new RealBigFloat(F)) {}
  Real(const char* s, const extLong& prec=defInputDigits) : RCReal(NULL) {
    constructFromString(s, prec);
  }
  Real(const std::string& s, const extLong& prec=defInputDigits) : RCReal(NULL){
    constructFromString(s.c_str(), prec);
  }

  /// \name Copy-Assignment-Destructor
  //@{
  /// copy constructor
  Real(const Real& rhs) : RCReal(rhs) {
    rep->incRef();
  }
  /// assignment operator
  Real& operator=(const Real& rhs) {
    if (this != &rhs) {
      rep->decRef();
      rep = rhs.rep;
      rep->incRef();
    }
    return *this;
  }
  /// destructor
  ~Real() {
    rep->decRef();
  }
  //@}

  /// \name Compound Assignment Operators
  //@{
  /// operator+=
  Real& operator+=(const Real& x);
  /// operator-=
  Real& operator-=(const Real& x);
  /// operator*=
  Real& operator*=(const Real& x);
  /// operator/=
  Real& operator/=(const Real& x);
  //@}

  /// \name Unary Minus, Increment and Decrement Operators
  //@{
  /// unary plus
  Real operator+() const {
    return Real(*this);
  }
  /// unary minus
  Real operator-() const {
    return -(*rep);
  }
  /// left increment operator (++i)
  Real& operator++() {
    *this += 1;
    return *this;
  }
  /// left decrement operator (--i)
  Real& operator--() {
    *this -= 1;
    return *this;
  }
  /// right increment operator (i++)
  Real operator++(int) {
    Real t(*this);
    *this += 1;
    return t;
  }
  /// right deccrement operator (i--)
  Real operator--(int) {
    Real t(*this);
    *this -= 1;
    return t;
  }
  //@}

  /// \name String Conversion Functions
  //@{
  /// set value from <tt>const char*</tt>
  void fromString(const char* s, const extLong& prec = defInputDigits) {
    *this = Real(s, prec);
  }
  /// convert to <tt>std::string</tt>
  /** give decimal string representation */
  std::string toString(long prec=defOutputDigits, bool sci=false) const {
    return rep->toString(prec, sci);
  }
  //@}

  /// \name Conversion Functions
  //@{
  /// convert to \c int
  int intValue() const {
    return static_cast<int>(longValue());
  }
  /// convert to \c long
  long longValue() const {
    return rep->longValue();
  }
  /// convert to \c float
  float floatValue() const {
    return static_cast<float>(doubleValue());
  }
  /// convert to \c double
  double doubleValue() const {
    return rep->doubleValue();
  }
  /// convert to \c BigInt
  BigInt BigIntValue() const {
    return rep->BigIntValue();
  }
  /// convert to \c BigRat
  BigRat BigRatValue() const {
    return rep->BigRatValue();
  }
  /// convert to \c BigFloat (approximate it first!)
  BigFloat BigFloatValue() const {
    return rep->BigFloatValue();
  }
  //@}

  /// \name Aprroximation Function
  //@{
  /// approximation
  Real approx(const extLong& r=defRelPrec, const extLong& a=defAbsPrec) const {
    return rep->approx(r, a);
  }
  //@}

  /// \name Helper Functions
  //@{
  /// sign function
  int sign() const {
    return rep->sgn();
  }
  /// isZero function
  bool isZero() const {
    return sign() == 0;
  }
  /// return true if interval contains zero
  bool isZeroIn() const {
    return rep->isZeroIn();
  }
  /// absolute value function
  Real abs() const {
    return (sign() >= 0) ? +(*this) : -(*this);
  }

  /// get mantissa of current approximate value
  BigInt getMantissa() const {
    return BigFloatValue().m();
  }
  /// get exponent of current approximate value
  long getExponent() const {
    return BigFloatValue().exp();
  }

  /// return true if error free otherwise return false;
  bool  isExact() const {
    return rep->isExact();
  }

  /// low bound of MSB
  extLong lMSB() const {
    return isExact() ? MSB():(rep->BigFloatValue()).lMSB();
  }
  /// upper bound of MSB
  extLong uMSB() const {
    return isExact() ? MSB():(rep->BigFloatValue()).uMSB();
  }
  /// MSB - Most Significant Bit
  extLong MSB() const {
    return rep->mostSignificantBit;
  }

  /// floor of log_2 of Error
  extLong flrLgErr() const {
    return rep->flrLgErr();
  }
  /// ceil of log_2 of Error
  extLong clLgErr() const {
    return rep->clLgErr();
  }

  /// division with desired precision
  Real div(const Real& x, const extLong& r) const;
  /// squareroot
  Real sqrt(const extLong& x) const {
    return rep->sqrt(x);
  }
  /// squareroot with initial approximation
  Real sqrt(const extLong& x, const BigFloat& A) const {
    return rep->sqrt(x, A);
  }

  /// correspond to the variables "u25, l25, v2p, v2m, v5p, v5m" in Expr
  void ULV_E(extLong &up, extLong &lp, extLong &v2p, extLong &v2m,
             extLong &v5p, extLong &v5m) const {
    rep->ULV_E(up, lp, v2p, v2m, v5p, v5m);
  }

  /// degree of polynomial P(x)
  unsigned long degree() const {
    return rep->degree();
  }
  /// \f$ lg(|| P(X) ||_2) \f$
  unsigned long length() const {
    return rep->length();
  }
  /// \f$ lg(|| P(X) ||_\infty) \f$
  unsigned long height() const {
    return rep->height();
  }
  //@}

  /// return Real(0)
  CGAL_CORE_EXPORT static const Real& getZero();
private:
  CGAL_CORE_EXPORT void constructFromString(const char *str, const extLong& prec);
};

#define CORE_REAL_ZERO Real::getZero()

const long halfLongMax = LONG_MAX /2;
const long halfLongMin = LONG_MIN /2;

struct _real_add {
  template <class T>
  static Real eval(const T& a, const T& b) {
    return a+b;
  }
  // specialized for two long values
  static Real eval(long a, long b) {
    if ((a > halfLongMax && b > halfLongMax) || (a < halfLongMin && b < halfLongMin))
      return BigInt(a)+BigInt(b);
    else
      return a+b;
  }
};

struct _real_sub {
  template <class T>
  static Real eval(const T& a, const T& b) {
    return a-b;
  }
  // specialized for two long values
  static Real eval(long a, long b) {
    if ((a > halfLongMax && b < halfLongMin) || (a < halfLongMin && b > halfLongMax))
      return BigInt(a)-BigInt(b);
    else
      return a-b;
  }
};

struct _real_mul {
  template <class T>
  static Real eval(const T& a, const T& b) {
    return a*b;
  }
  // specialized for two long values
  static Real eval(long a, long b) {
    if (flrLg(a) + flrLg(b) >= static_cast<int>(LONG_BIT-2))
      return BigInt(a)*BigInt(b);
    else
      return a*b;
  }
};

template <class Op>
struct _real_binary_op {
  static Real eval(const RealRep& a, const RealRep& b) {
    if (a.ID() == REAL_BIGRAT || b.ID() == REAL_BIGRAT) {
      if (!a.isExact()) { // a must be a BigFloat and b must be a BigRat
        BigFloat bf_a = a.BigFloatValue(), bf_b;
        bf_b.approx(b.BigRatValue(), CORE_posInfty, -bf_a.flrLgErr());
        return Op::eval(bf_a, bf_b);
      } else if (!b.isExact()) { // a must be a BigRat and b must be a BigFloat
        BigFloat bf_a, bf_b = b.BigFloatValue();
        bf_a.approx(a.BigRatValue(), CORE_posInfty, -bf_b.flrLgErr());
        return Op::eval(bf_a, bf_b);
      } else // both are BigRat
        return Op::eval(a.BigRatValue(), b.BigRatValue());
    } else if (a.ID() == REAL_BIGFLOAT || b.ID() == REAL_BIGFLOAT
               || a.ID() == REAL_DOUBLE || b.ID() == REAL_DOUBLE) {
      return Op::eval(a.BigFloatValue(), b.BigFloatValue());
    } else if (a.ID() == REAL_BIGINT || b.ID() == REAL_BIGINT) {
      return Op::eval(a.BigIntValue(), b.BigIntValue());
    } else { // a.ID() == REAL_LONG && b.ID() == REAL_LONG
      return Op::eval(a.longValue(), b.longValue());
    }
  }
};

typedef _real_binary_op<_real_add> real_add;
typedef _real_binary_op<_real_sub> real_sub;
typedef _real_binary_op<_real_mul> real_mul;

struct real_div {
  static Real eval(const RealRep& a, const RealRep& b, const extLong& r) {
    if (a.ID() == REAL_BIGRAT || b.ID() == REAL_BIGRAT) {
      if (!a.isExact()) { // a must be a BigFloat and b must be a BigRat
        BigFloat bf_a = a.BigFloatValue(), bf_b;
        bf_b.approx(b.BigRatValue(), bf_a.MSB() - bf_a.flrLgErr() + 1, CORE_posInfty);
        return bf_a.div(bf_b, r);
      } else if (!b.isExact()) { // a must be a BigRat and b must be a BigFloat
        BigFloat bf_a, bf_b = b.BigFloatValue();
        bf_a.approx(a.BigRatValue(), bf_b.MSB() - bf_b.flrLgErr() + 1, CORE_posInfty);
        return bf_a.div(bf_b, r);
      } else // both are BigRat
        return a.BigRatValue()/b.BigRatValue();
    } else if (a.ID() == REAL_BIGFLOAT || b.ID() == REAL_BIGFLOAT
               || a.ID() == REAL_DOUBLE || b.ID() == REAL_DOUBLE) {
      return a.BigFloatValue().div(b.BigFloatValue(), r);
    } else if (a.ID() == REAL_BIGINT || b.ID() == REAL_BIGINT) {
      return BigRat(a.BigIntValue(), b.BigIntValue());
    } else { // a.ID() == REAL_LONG && b.ID() == REAL_LONG
      return BigRat(a.longValue(), b.longValue());
    }
  }
};

CGAL_CORE_EXPORT std::istream& operator>>(std::istream& i, Real& r);

inline std::ostream& operator<<(std::ostream& o, const Real& r) {
  return r.getRep().operator<<(o);
}

inline Real& Real::operator+=(const Real& rhs) {
  *this = real_add::eval(getRep(), rhs.getRep());
  return *this;
}
inline Real& Real::operator-=(const Real& rhs) {
  *this = real_sub::eval(getRep(), rhs.getRep());
  return *this;
}
inline Real& Real::operator*=(const Real& rhs) {
  *this = real_mul::eval(getRep(), rhs.getRep());
  return *this;
}
inline Real& Real::operator/=(const Real& rhs) {
  *this = real_div::eval(getRep(), rhs.getRep(), defRelPrec);
  return *this;
}

// operator+
inline Real operator+(const Real& x, const Real& y) {
  return real_add::eval(x.getRep(), y.getRep());
}
// operator-
inline Real operator-(const Real& x, const Real& y) {
  return real_sub::eval(x.getRep(), y.getRep());
}
// operator*
inline Real operator*(const Real& x, const Real& y) {
  return real_mul::eval(x.getRep(), y.getRep());
}
// operator/
inline Real operator/(const Real& x, const Real& y) {
  return real_div::eval(x.getRep(), y.getRep(), defRelPrec);
}
// div w/ precision
inline Real Real::div(const Real& x, const extLong& r) const {
  return real_div::eval(getRep(), x.getRep(), r);
}

inline int cmp(const Real& x, const Real& y) {
  return (x-y).sign();
}
inline bool operator==(const Real& x, const Real& y) {
  return cmp(x, y) == 0;
}
inline bool operator!=(const Real& x, const Real& y) {
  return cmp(x, y) != 0;
}
inline bool operator>=(const Real& x, const Real& y) {
  return cmp(x, y) >= 0;
}
inline bool operator>(const Real& x, const Real& y) {
  return cmp(x, y) > 0;
}
inline bool operator<=(const Real& x, const Real& y) {
  return cmp(x, y) <= 0;
}
inline bool operator<(const Real& x, const Real& y) {
  return cmp(x, y) < 0;
}

/// floor function
CGAL_CORE_EXPORT BigInt floor(const Real&, Real&);
/// power function
CGAL_CORE_EXPORT Real pow(const Real&, unsigned long);

/// return sign
inline int sign(const Real& r) {
  return r.sign();
}
/// is zero?
inline bool isZero(const Real& r) {
  return r.sign() == 0;
}
/// absolute value
inline Real abs(const Real& x) {
  return x.abs();
}
/// absolute value (same as abs)
inline Real fabs(const Real& x) {
  return abs(x);
}
/// floor
inline BigInt floor(const Real& r) {
  Real tmp;
  return floor(r, tmp);
}
/// ceiling
inline BigInt ceil(const Real& r) {
  return -floor(-r);
}
/// power
inline Real power(const Real& r, unsigned long p) {
  return pow(r, p);
}
/// square root
inline Real sqrt(const Real& x) {
  return x.sqrt(defAbsPrec);
}

// class Realbase_for (need defined after Real)
// unary minus operator
template <class T>
inline Real Realbase_for<T>::operator-() const {
  return -ker;
}
template <>
inline Real RealLong::operator-() const {
  return ker < -LONG_MAX ? -BigInt(ker) : -ker;
}

} //namespace CORE
#endif // _CORE_REAL_H_
