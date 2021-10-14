/****************************************************************************
 * Core Library Version 1.7, August 2004
 * Copyright (c) 1995-2004 Exact Computation Project
 * All rights reserved.
 *
 * This file is part of CGAL (www.cgal.org).
 *
 * File: BigRat.h
 * Synopsis:
 *                 a wrapper class for mpq from GMP
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

#ifndef _CORE_BIGRAT_H_
#define _CORE_BIGRAT_H_

#include <CGAL/CORE/BigInt.h>

namespace CORE {

#ifdef CGAL_CORE_USE_GMP_BACKEND
  typedef boost::multiprecision::mpq_rational Q;
#else
  typedef  boost::multiprecision::cpp_rational Q;
#endif

class BigRatRep : public RCRepImpl<BigRatRep> {
public:
  BigRatRep()
  : mp()
  {}

  // Note : should the copy-ctor be alloed at all ? [Sylvain Pion]
  BigRatRep(const BigRatRep& z)  : RCRepImpl<BigRatRep>(), mp(z.mp)
  {}


  BigRatRep(signed char c)
    : mp(c)
  {}

  BigRatRep(unsigned char c)
    : mp(c)
  {}

  BigRatRep(signed int i)
    : mp(i)
  {}

  BigRatRep(unsigned int i)
    : mp(i)
  {}

  BigRatRep(signed short int s)
    : mp(s)
  {}

  BigRatRep(unsigned short int s)
    : mp(s)
  {}

  BigRatRep(signed long int l)
    : mp(l)
  {}

  BigRatRep(unsigned long int l)
    : mp(l)
  {}

  BigRatRep(float f)
    : mp(f)
  {}

  BigRatRep(double d)
    : mp(d)
  {}

  BigRatRep(const char* s)
    : mp(s)
  {}

  BigRatRep(const std::string& s)
    : mp(s)
  {}

  explicit BigRatRep(const Z& q)
  : mp(q)
  {}

  BigRatRep(const Z& n, const Z& d)
    : mp(n,d)
  {}


  //CGAL_CORE_EXPORT CORE_NEW(BigRatRep)
//  CGAL_CORE_EXPORT CORE_DELETE(BigRatRep)

  const Q& get_mp() const {
    return mp;
  }
  Q& get_mp() {
    return mp;
  }
private:
  Q mp;
}; //BigRatRep

class BigFloat;

typedef RCImpl<BigRatRep> RCBigRat;
class BigRat : public RCBigRat {
public:
  /// \name Constructors
  //@{
  /// default constructor
  BigRat() : RCBigRat(new BigRatRep()) {}
  /// constructor for <tt>signed char</tt>
  BigRat(signed char x) : RCBigRat(new BigRatRep(x)) {}
  /// constructor for <tt>unsigned char</tt>
  BigRat(unsigned char x) : RCBigRat(new BigRatRep(x)) {}
  /// constructor for <tt>signed short int</tt>
  BigRat(signed short int x) : RCBigRat(new BigRatRep(x)) {}
  /// constructor for <tt>unsigned short int</tt>
  BigRat(unsigned short int x) : RCBigRat(new BigRatRep(x)) {}
  /// constructor for <tt>signed int</tt>
  BigRat(signed int x) : RCBigRat(new BigRatRep(x)) {}
  /// constructor for <tt>unsigned int</tt>
  BigRat(unsigned int x) : RCBigRat(new BigRatRep(x)) {}
  /// constructor for <tt>signed long int</tt>
  BigRat(signed long int x) : RCBigRat(new BigRatRep(x)) {}
  /// constructor for <tt>unsigned long int</tt>
  BigRat(unsigned long int x) : RCBigRat(new BigRatRep(x)) {}
  /// constructor for <tt>float</tt>
  BigRat(float x) : RCBigRat(new BigRatRep(x)) {}
  /// constructor for <tt>double</tt>
  BigRat(double x) : RCBigRat(new BigRatRep(x)) {}
  /// constructor for <tt>const char*</tt> with base
  BigRat(const char* s) : RCBigRat(new BigRatRep(s)) {}
  /// constructor for <tt>std::string</tt> with base
  BigRat(const std::string& s) : RCBigRat(new BigRatRep(s)) {}
  /// constructor for <tt>mpq_srcptr</tt>
  explicit BigRat(const Z& z) : RCBigRat(new BigRatRep(z)) {}
  /// constructor for <tt>BigInt</tt>
  BigRat(const BigInt& z) : RCBigRat(new BigRatRep(z.get_mp())) {}
  /// constructor for two <tt>BigInts</tt>
  BigRat(const BigInt& n, const BigInt& d)
      : RCBigRat(new BigRatRep(n.get_mp(), d.get_mp())) {}
  /// constructor for <tt>BigFloat</tt>
  BigRat(const BigFloat&);
  //@}

  /// \name Copy-Assignment-Destructor
  //@{
  /// copy constructor
  BigRat(const BigRat& rhs) : RCBigRat(rhs) {
    rep->incRef();
  }
  /// assignment operator
  BigRat& operator=(const BigRat& rhs) {
    if (this != &rhs) {
      rep->decRef();
      rep = rhs.rep;
      rep->incRef();
    }
    return *this;
  }
  /// destructor
  ~BigRat() {
    rep->decRef();
  }
  //@}

  /// \name Overloaded operators
  //@{
  BigRat& operator +=(const BigRat& rhs) {
    makeCopy();
    get_mp() += rhs.get_mp();
    return *this;
  }
  BigRat& operator -=(const BigRat& rhs) {
    makeCopy();
    get_mp() -= rhs.get_mp();
    return *this;
  }
  BigRat& operator *=(const BigRat& rhs) {
    makeCopy();
    get_mp() *= rhs.get_mp();
    return *this;
  }
  BigRat& operator /=(const BigRat& rhs) {
    makeCopy();
    get_mp() /= rhs.get_mp();
    return *this;
  }
  BigRat& operator <<=(unsigned long ul) {
    makeCopy();
    assert(false);
    // AF no shift for Q   get_mp() <<= ul;
    return *this;
  }
  BigRat& operator >>=(unsigned long ul) {
    makeCopy();
    assert(false);
    // AF no >> for Q    get_mp() >>= ul;
    return *this;
  }
  //@}

  /// \name div2, unary, increment, decrement operators
  //@{

  /// exact division by 2 (this method is provided for compatibility)
  BigRat div2() const {
    BigRat r; BigRat t(2);     // probably not most efficient way
    r.get_mp() = get_mp() / t.get_mp();
    return r;
  }
  BigRat operator+() const {
    return BigRat(*this);
  }
  BigRat operator-() const {
    BigRat r;
    r.get_mp() = - get_mp();
    return r;
  }
  BigRat& operator++() {
    makeCopy();
    ++get_mp();
    return *this;
  }
  BigRat& operator--() {
    makeCopy();
    --get_mp();
    return *this;
  }
  BigRat operator++(int) {
    BigRat r(*this);
    ++(*this);
    return r;
  }
  BigRat operator--(int) {
    BigRat r(*this);
    --(*this);
    return r;
  }
  //@}

  /// \name Helper Functions
  //@{
  /// Canonicalize
  void canonicalize() {
    makeCopy();
    assert(false);
    // AF todo    mpq_canonicalize(get_mp());
  }
  /// Has Exact Division
  static bool hasExactDivision() {
    return true;
  }

  /// return mpz pointer of numerator (const)
  Z get_num_mp() const {
    return numerator(get_mp());
  }
  /// return mpz pointer of numerator  // no references as numerator() returns a copy
  Z get_num_mp() {
    return numerator(get_mp());
  }
  /// return mpz pointer of denominator
  Z get_den_mp() const {
    return denominator(get_mp());
  }
  /// return mpz pointer of denominator
  Z get_den_mp() {
    return denominator(get_mp());
  }

  /// get mpq pointer (const)
  const Q& get_mp() const {
    return rep->get_mp();
  }
  /// get mpq pointer
  Q& get_mp() {
    return rep->get_mp();
  }
  //@}

  /// \name String Conversion Functions
  //@{
  /// set value from <tt>const char*</tt>
  int set_str(const char* s, int base = 0) {
    makeCopy();
    get_mp() = Q(s);
    return 0;
  }
  /// convert to <tt>std::string</tt>
  std::string get_str(int base = 10) const {
    return get_mp().convert_to<std::string>();
  }
  //@}

  /// \name Conversion Functions
  //@{
  /// intValue
  int intValue() const {
    return get_mp().convert_to<int>();
  }
  /// longValue
  long longValue() const {
    return get_mp().convert_to<long>();
  }

  /// doubleValue
  double doubleValue() const {
    return get_mp().convert_to<double>();
  }
  /// BigIntValue
  BigInt BigIntValue() const {
    BigInt r;
    Z rem;
    divide_qr(get_num_mp(),get_den_mp(), r.get_mp(), rem);
    return r;
  }
  //@}
}; //BigRat class

inline BigRat operator+(const BigRat& a, const BigRat& b) {
  BigRat r;
  r.get_mp() = a.get_mp() + b.get_mp();
  return r;
}
inline BigRat operator-(const BigRat& a, const BigRat& b) {
  BigRat r;
  r.get_mp() = a.get_mp() - b.get_mp();
  return r;
}
inline BigRat operator*(const BigRat& a, const BigRat& b) {
  BigRat r;
  r.get_mp() = a.get_mp() * b.get_mp();
  return r;
}
inline BigRat operator/(const BigRat& a, const BigRat& b) {
  BigRat r;
  r.get_mp() = a.get_mp() / b.get_mp();
  return r;
}
// Chee (3/19/2004):
//   The following definitions of div_exact(x,y) and gcd(x,y)
//   ensures that in Polynomial<NT>
/// divisible(x,y) = "x | y"
inline BigRat div_exact(const BigRat& x, const BigRat& y) {
        BigRat z;
        z.get_mp() =  x.get_mp() / y.get_mp();
        return z;
}
/// numerator
inline BigInt numerator(const BigRat& a) {
  return BigInt(a.get_num_mp());
}
/// denominator
inline BigInt denominator(const BigRat& a) {
  return BigInt(a.get_den_mp());
}

inline BigRat gcd(const BigRat& x, const BigRat& y) {
  //        return BigRat(1);  // Remark: we may want replace this by
                           // the definition of gcd of a quotient field
                           // of a UFD [Yap's book, Chap.3]
  //Here is one possible definition: gcd of x and y is just the
  //gcd of the numerators of x and y divided by the gcd of the
  //denominators of x and y.
  BigInt n = gcd(numerator(x), numerator(y));
  BigInt d = gcd(denominator(x), denominator(y));
  return BigRat(n,d);

}
// Chee: 8/8/2004: need isDivisible to compile Polynomial<BigRat>
// A trivial implementation is to return true always. But this
// caused tPolyRat to fail.
// So we follow the definition of
// Expr::isDivisible(e1, e2) which checks if e1/e2 is an integer.
inline bool isInteger(const BigRat& x) {
  return BigInt(x.get_den_mp()) == 1;
}
inline bool isDivisible(const BigRat& x, const BigRat& y) {
  BigRat r;
  r.get_mp() = x.get_mp() / y.get_mp();
  return isInteger(r);
}
inline BigRat operator<<(const BigRat& a, unsigned long ul) {
  BigRat r;
  assert(false);
  //AF todo   no << for Q     r.get_mp() = a.get_mp() << int(ul) ;
  return r;
}
inline BigRat operator>>(const BigRat& a, unsigned long ul) {
  BigRat r;
  assert(false);
  // AF todo no >> for Q    r.get_mp() =  a.get_mp() >> ul;
  return r;
}

inline int cmp(const BigRat& x, const BigRat& y) {
  return x.get_mp().compare(y.get_mp());
}
inline bool operator==(const BigRat& a, const BigRat& b) {
  return cmp(a, b) == 0;
}
inline bool operator!=(const BigRat& a, const BigRat& b) {
  return cmp(a, b) != 0;
}
inline bool operator>=(const BigRat& a, const BigRat& b) {
  return cmp(a, b) >= 0;
}
inline bool operator>(const BigRat& a, const BigRat& b) {
  return cmp(a, b) > 0;
}
inline bool operator<=(const BigRat& a, const BigRat& b) {
  return cmp(a, b) <= 0;
}
inline bool operator<(const BigRat& a, const BigRat& b) {
  return cmp(a, b) < 0;
}

inline std::ostream& operator<<(std::ostream& o, const BigRat& x) {
  return o << x.get_mp();
}
inline std::istream& operator>>(std::istream& i, BigRat& x) {
  x.makeCopy();
  return i >> x.get_mp();
}

/// sign
inline int sign(const BigRat& a) {
  return sign(a.get_mp());
}
/// abs
inline BigRat abs(const BigRat& a) {
  BigRat r;
  r.get_mp() = abs( a.get_mp());
  return r;
}
/// neg
inline BigRat neg(const BigRat& a) {
  BigRat r;
  r.get_mp() = - a.get_mp();
  return r;
}
/// div2
inline BigRat div2(const BigRat& a) {
  BigRat r(a);
  return r.div2();
}
/// longValue
inline long longValue(const BigRat& a) {
  return a.longValue();
}
/// doubleValue
inline double doubleValue(const BigRat& a) {
  return a.doubleValue();
}
/// return BigInt value
inline BigInt BigIntValue(const BigRat& a) {
  return a.BigIntValue();
}


} //namespace CORE
#endif // _CORE_BIGRAT_H_
