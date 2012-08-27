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
 * File: RealRep.h
 * Synopsis: 
 * 		Internal Representation for Real
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
#ifndef _CORE_REALREP_H_
#define _CORE_REALREP_H_
#include "BigFloat.h"

namespace CORE { 

class Real;

class RealRep {
public:
  extLong mostSignificantBit;
public:
  RealRep() : refCount(1) {}
  virtual ~RealRep() {}
  virtual int ID() const = 0;

  virtual long longValue() const = 0;
  virtual double doubleValue() const = 0;
  virtual BigInt BigIntValue() const = 0;
  virtual BigRat BigRatValue() const = 0;
  virtual BigFloat BigFloatValue() const = 0;

  virtual BigFloat approx(const extLong&, const extLong&) const = 0;
  virtual Real operator-() const = 0;

  virtual bool isExact() const = 0;
  virtual int sgn() const = 0;
  virtual bool isZeroIn() const = 0;

  virtual BigFloat sqrt(const extLong&) const = 0;
  virtual BigFloat sqrt(const extLong&, const BigFloat&) const = 0;

  virtual void ULV_E(extLong &, extLong&, extLong&,
		  extLong&, extLong&, extLong&) const = 0;
  virtual extLong flrLgErr() const = 0;
  virtual extLong clLgErr() const = 0;
  virtual unsigned long degree() const = 0;
  virtual unsigned long length() const = 0;
  virtual unsigned long height() const = 0;

  virtual std::string toString(long prec, bool sci) const = 0;
  virtual std::ostream& operator<<(std::ostream& o) const = 0;
public:
  void incRef() {
    ++refCount;
  }
  void decRef() {
    if (--refCount == 0)
      delete this;
  }
  int getRefCount() const {
    return refCount;
  }
private:
  int refCount;
};//realRep class

template <class T>
class Realbase_for : public RealRep {
public:
  CORE_MEMORY(Realbase_for)
  Realbase_for(const T& k);
  ~Realbase_for() {}
  int ID() const;

  long longValue() const {
    return ker.longValue();
  }
  double doubleValue() const {
    return ker.doubleValue();
  }
  BigInt BigIntValue() const {
    return BigInt(ker);
  }
  BigRat BigRatValue() const {
    return BigRat(ker);
  }
  BigFloat BigFloatValue() const {
    return BigFloat(ker);
  }

  BigFloat approx(const extLong&, const extLong&) const;
  Real operator-() const;

  bool isExact() const {
    return true;
  }
  int sgn() const {
    return ker > 0.0 ? 1 : ( ker == 0.0 ? 0 : -1);
  }
  bool isZeroIn() const {
    return ker == 0.0;
  }

  BigFloat sqrt(const extLong&) const;
  BigFloat sqrt(const extLong&, const BigFloat&) const;

  void ULV_E(extLong &, extLong&, extLong&, extLong&, extLong&, extLong&) const;
  extLong flrLgErr() const {
    return CORE_negInfty;
  }
  extLong clLgErr() const {
    return CORE_negInfty;
  }
  unsigned long degree() const {
    return 1;
  }
  unsigned long length() const;
  unsigned long height() const;

  std::string toString(long, bool) const {
    std::stringstream st;
    st << ker;
    return st.str();
  }
  std::ostream& operator<<(std::ostream& o) const {
    return o << ker;
  }
private:
  T ker;
};//Realbase_for class

typedef Realbase_for<long> RealLong;
typedef Realbase_for<double> RealDouble;
typedef Realbase_for<BigInt> RealBigInt;
typedef Realbase_for<BigRat> RealBigRat;
typedef Realbase_for<BigFloat> RealBigFloat;

enum { REAL_LONG, REAL_DOUBLE, REAL_BIGINT, REAL_BIGRAT, REAL_BIGFLOAT };

// constructors
template<>
inline RealLong::Realbase_for(const long& l) : ker(l) {
  mostSignificantBit = (ker != 0 ) ? extLong(flrLg(ker)) : CORE_negInfty;
}
template<>
inline RealDouble::Realbase_for(const double& d) : ker(d) {
  mostSignificantBit = BigFloat(ker).MSB();
}
template<>
inline RealBigInt::Realbase_for(const BigInt& l) : ker(l) {
  mostSignificantBit = (sign(ker)) ? extLong(floorLg(ker)) : CORE_negInfty;
}
template<>
inline RealBigRat::Realbase_for(const BigRat& l) : ker(l) {
  mostSignificantBit = BigFloat(ker).MSB();
}
template<>
inline RealBigFloat::Realbase_for(const BigFloat& l) : ker(l) {
  mostSignificantBit = ker.MSB();
}

// ID()
template<>
inline int RealLong::ID() const {
  return REAL_LONG;
}
template<>
inline int RealDouble::ID() const {
  return REAL_DOUBLE;
}
template<>
inline int RealBigInt::ID() const {
  return REAL_BIGINT;
}
template<>
inline int RealBigRat::ID() const {
  return REAL_BIGRAT;
}
template<>
inline int RealBigFloat::ID() const {
  return REAL_BIGFLOAT;
}

// cast functions
template<>
inline long RealLong::longValue() const {
  return ker;
}
template<>
inline long RealDouble::longValue() const {
  return static_cast<long>(ker);
}
template<>
inline double RealLong::doubleValue() const {
  return static_cast<double>(ker);
}
template<>
inline double RealDouble::doubleValue() const {
  return ker;
}
template<>
inline BigInt   RealBigInt::BigIntValue() const {
  return ker;
}
template<>
inline BigInt   RealBigRat::BigIntValue() const {
  return ker.BigIntValue();
}
template<>
inline BigInt RealBigFloat::BigIntValue() const {
  return ker.BigIntValue();
}
template<>
inline BigRat   RealBigRat::BigRatValue() const {
  return ker;
}
template<>
inline BigRat RealBigFloat::BigRatValue() const {
  return ker.BigRatValue();
}
template<>
inline BigFloat RealBigFloat::BigFloatValue() const {
  return ker;
}

// isExact()
template<>
inline bool RealBigFloat::isExact() const {
  return ker.isExact();
}

// sign()
template<>
inline int RealBigInt::sgn() const {
  return sign(ker);
}
template<>
inline int RealBigRat::sgn() const {
  return sign(ker);
}
template<>
inline int RealBigFloat::sgn() const {
  return ker.sign();
}

// isZeroIn()
template<>
inline bool RealBigInt::isZeroIn() const {
  return sign(ker) == 0;
}
template<>
inline bool RealBigRat::isZeroIn() const {
  return sign(ker) == 0;
}
template<>
inline bool RealBigFloat::isZeroIn() const {
  return ker.isZeroIn();
}

// approx
template <class T>
inline BigFloat Realbase_for<T>::approx(const extLong& r, const extLong& a) const {
  BigFloat x;
  x.approx(ker, r, a);
  return x;
}
template <>
inline BigFloat RealLong::approx(const extLong& r, const extLong& a) const {
  BigFloat x;
  x.approx(BigInt(ker), r, a);
  return x;
}
template <>
inline BigFloat RealDouble::approx(const extLong& r, const extLong& a) const {
  BigFloat x;
  x.approx(BigRat(ker), r, a);
  return x;
}

// sqrt
template <class T>
inline BigFloat Realbase_for<T>::sqrt(const extLong& a) const {
  return BigFloat(ker).sqrt(a);
}
template <class T>
inline BigFloat Realbase_for<T>::sqrt(const extLong& a, const BigFloat& A) const {
  return BigFloat(ker).sqrt(a, A);
}

// ULV_E()
template<>
inline void RealLong::ULV_E(extLong &up, extLong &lp, extLong &v2p,
                            extLong &v2m, extLong &v5p, extLong &v5m) const {
  // TODO : extract the power of 5.
  up = lp = v2p = v2m = v5p = v5m = EXTLONG_ZERO;
  if (ker == 0)
    return;

  // Extract the power of 2.
  unsigned long exp = 0;
  unsigned long tmp_ker = ker;
  while ((tmp_ker&1) != 0) {
    tmp_ker = tmp_ker/2;
    ++exp;
  }
  up = clLg(tmp_ker);
  lp = 0;
  v2p = exp;
}
template<>
inline void RealDouble::ULV_E(extLong &up, extLong &lp, extLong &v2p,
                              extLong &v2m, extLong &v5p, extLong &v5m) const {
  // TODO : can probably be made faster using frexp() or such.
  // TODO : extract the power of 5.
  BigRat R = BigRat(ker);
  up  = ceilLg(numerator(R));
  v2m = ceilLg(denominator(R));
  lp = v2p = v5m = v5p = EXTLONG_ZERO;
}
template<>
inline void RealBigInt::ULV_E(extLong &up, extLong &lp, extLong &v2p,
                              extLong &v2m, extLong &v5p, extLong &v5m) const {
  up = lp = v2p = v2m = v5p = v5m = EXTLONG_ZERO;
  if (ker == 0)
    return;

  // Extract power of 5.
  int exp5;
  BigInt remainder5;
  getKaryExpo(ker, remainder5, exp5, 5);
  v5p = exp5;
  // Extract power of 2.
  int exp2 = getBinExpo(remainder5);
  up = ceilLg(remainder5) - exp2;
  v2p = exp2;
}
template<>
inline void RealBigRat::ULV_E(extLong &up, extLong &lp, extLong &v2p,
                              extLong &v2m, extLong &v5p, extLong &v5m) const {
  up = lp = v2p = v2m = v5p = v5m = EXTLONG_ZERO;
  if (ker == 0)
    return;

  // Extract power of 5.
  int exp5;
  BigInt num5, den5;
  getKaryExpo(numerator(ker), num5, exp5, 5);
  if (exp5 != 0) {
    v5p = exp5;
    den5 = denominator(ker);
  } else {
    getKaryExpo(denominator(ker), den5, exp5, 5);
    v5m = exp5;
  }

  // Now we work with num5/den5.
  int exp2 = getBinExpo(num5);
  if (exp2 != 0) {
    v2p = exp2;
  } else {
    exp2 = getBinExpo(den5);
    v2m = exp2;
  }

  up = ceilLg(num5) - v2p;
  lp = ceilLg(den5) - v2m;
}
template<>
inline void RealBigFloat::ULV_E(extLong &up, extLong &lp, extLong &v2p,
                                extLong &v2m, extLong &v5p, extLong &v5m) const {
  // TODO : extract power of 5.
  up = lp = v2p = v2m = v5p = v5m = EXTLONG_ZERO;
  BigRat R = ker.BigRatValue();
  up  = ceilLg(numerator(R));
  v2m = ceilLg(denominator(R));
}

// flrLgErr && clLgErr
template<>
inline extLong RealBigFloat::flrLgErr() const {
  return ker.flrLgErr();
}
template<>
inline extLong RealBigFloat::clLgErr() const {
  return ker.clLgErr();
}

// height && length
template<>
inline unsigned long RealLong::length() const {
  return clLg(1+ core_abs(ker));
}	// length is (log_2(1+ker^2)) /2.

template<>
inline unsigned long RealLong::height() const {
  return clLg(core_max(1L, core_abs(ker)));
}	// height is max{1, |ker|}

template<>
inline unsigned long RealDouble::length() const {
  BigRat R  = BigRat(ker);
  long ln = 1 + ceilLg(numerator(R));
  long ld = 1 + ceilLg(denominator(R));
  return (ln>ld) ? ln : ld; ///< an upper bound on log_2(sqrt(num^2+den^2))
}

template<>
inline unsigned long RealDouble::height() const {
  BigRat R  = BigRat(ker);
  long ln = ceilLg(numerator(R));
  long ld = ceilLg(denominator(R));
  return (ln>ld) ? ln : ld; ///< an upper bound on log_2(max(|num|, |den|))
}
template<>
inline unsigned long RealBigInt::length() const {
  return ceilLg(1 + abs(ker));
}

template<>
inline unsigned long RealBigInt::height() const {
  BigInt r(abs(ker));
  if (r<1)
    r = 1;
  return ceilLg(r);
}

template<>
inline unsigned long RealBigFloat::length() const {
  // Chen Li: A bug fixed.
  // The statement in the older version with the bug was:
  //   BigRat R  = BigRat(ker);
  // The BigRat(BigFloat) actually is a
  // conversion operator (defined in BigFloat.h), _NOT_
  // an ordinary class constructor! The C++ language
  // specify that an intialization is not an assignment
  // but a constructor operation!
  // Considering that BigRat(BigFloat) is a conversion
  // operator not really a constructor. The programmer's
  // intent is obvious to do an assignment.
  // However, the g++ seems to be confused by the above
  // initialization.
  BigRat R  = ker.BigRatValue();
  long   ln = 1 + ceilLg(numerator(R));
  long   ld = 1 + ceilLg(denominator(R));
  return ( ln > ld ) ? ln : ld;
}

template<>
inline unsigned long RealBigFloat::height() const {
  // Chen Li: A bug fixed. The old statement with the bug was:
  //   BigRat R  = BigRat(ker);
  // Detailed reasons see above (in RealBigFloat::length()!
  BigRat R  = ker.BigRatValue();
  long     ln = ceilLg(numerator(R));
  long     ld = ceilLg(denominator(R));
  return   ( ln > ld ) ? ln : ld;
}

template<>
inline unsigned long RealBigRat::length() const {
  long ln = 1 + ceilLg(numerator(ker));
  long ld = 1 + ceilLg(denominator(ker));
  return ( ln > ld ) ? ln : ld;
}

template<>
inline unsigned long RealBigRat::height() const {
  long ln = ceilLg(numerator(ker));
  long ld = ceilLg(denominator(ker));
  return (ln > ld ) ? ln : ld;
}

// toString()
template<>
inline std::string RealBigInt::toString(long, bool) const {
  return ker.get_str();
}
template<>
inline std::string RealBigRat::toString(long, bool) const {
  return ker.get_str();
}
template<>
inline std::string RealBigFloat::toString(long prec, bool sci) const {
  return ker.toString(prec, sci);
}

} //namespace CORE
#endif // _CORE_REALREP_H_
