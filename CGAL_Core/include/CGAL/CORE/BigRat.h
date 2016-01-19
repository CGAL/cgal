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
 * File: BigRat.h
 * Synopsis: 
 * 		a wrapper class for mpq from GMP
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
 ***************************************************************************/

#ifndef _CORE_BIGRAT_H_
#define _CORE_BIGRAT_H_

#include <CGAL/CORE/BigInt.h>

namespace CORE { 

class BigRatRep : public RCRepImpl<BigRatRep> {
public:
  BigRatRep() {
    mpq_init(mp);
  }
  // Note : should the copy-ctor be alloed at all ? [Sylvain Pion]
  BigRatRep(const BigRatRep& z)  : RCRepImpl<BigRatRep>() {
    mpq_init(mp);
    mpq_set(mp, z.mp);
  }
  BigRatRep(signed char c) {
    mpq_init(mp);
    mpq_set_si(mp, c, 1);
  }
  BigRatRep(unsigned char c) {
    mpq_init(mp);
    mpq_set_ui(mp, c, 1);
  }
  BigRatRep(signed int i) {
    mpq_init(mp);
    mpq_set_si(mp, i, 1);
  }
  BigRatRep(unsigned int i) {
    mpq_init(mp);
    mpq_set_ui(mp, i, 1);
  }
  BigRatRep(signed short int s) {
    mpq_init(mp);
    mpq_set_si(mp, s, 1);
  }
  BigRatRep(unsigned short int s) {
    mpq_init(mp);
    mpq_set_ui(mp, s, 1);
  }
  BigRatRep(signed long int l) {
    mpq_init(mp);
    mpq_set_si(mp, l, 1);
  }
  BigRatRep(unsigned long int l) {
    mpq_init(mp);
    mpq_set_ui(mp, l, 1);
  }
  BigRatRep(float f) {
    mpq_init(mp);
    mpq_set_d(mp, f);
  }
  BigRatRep(double d) {
    mpq_init(mp);
    mpq_set_d(mp, d);
  }
  BigRatRep(const char* s) {
    mpq_init(mp);
    mpq_set_str(mp, s, 0);
  }
  BigRatRep(const std::string& s) {
    mpq_init(mp);
    mpq_set_str(mp, s.c_str(), 0);
  }
  explicit BigRatRep(mpq_srcptr q) {
    mpq_init(mp);
    mpq_set(mp, q);
  }
  BigRatRep(mpz_srcptr z) {
    mpq_init(mp);
    mpq_set_z(mp, z);
  }
  BigRatRep(mpz_srcptr n, mpz_srcptr d) {
    mpq_init(mp);
    mpz_set(mpq_numref(mp), n);
    mpz_set(mpq_denref(mp), d);
    mpq_canonicalize(mp);
  }
  ~BigRatRep() {
    mpq_clear(mp);
  }

  CGAL_CORE_EXPORT CORE_NEW(BigRatRep)
  CGAL_CORE_EXPORT CORE_DELETE(BigRatRep)

  mpq_srcptr get_mp() const {
    return mp;
  }
  mpq_ptr get_mp() {
    return mp;
  }
private:
  mpq_t mp;
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
  explicit BigRat(mpq_srcptr z) : RCBigRat(new BigRatRep(z)) {}
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
    mpq_add(get_mp(), get_mp(), rhs.get_mp());
    return *this;
  }
  BigRat& operator -=(const BigRat& rhs) {
    makeCopy();
    mpq_sub(get_mp(), get_mp(), rhs.get_mp());
    return *this;
  }
  BigRat& operator *=(const BigRat& rhs) {
    makeCopy();
    mpq_mul(get_mp(), get_mp(), rhs.get_mp());
    return *this;
  }
  BigRat& operator /=(const BigRat& rhs) {
    makeCopy();
    mpq_div(get_mp(), get_mp(), rhs.get_mp());
    return *this;
  }
  BigRat& operator <<=(unsigned long ul) {
    makeCopy();
    mpq_mul_2exp(get_mp(), get_mp(), ul);
    return *this;
  }
  BigRat& operator >>=(unsigned long ul) {
    makeCopy();
    mpq_div_2exp(get_mp(), get_mp(), ul);
    return *this;
  }
  //@}

  /// \name div2, unary, increment, decrement operators
  //@{

  /// exact division by 2 (this method is provided for compatibility)
  BigRat div2() const {
    BigRat r; BigRat t(2);     // probably not most efficient way
    mpq_div(r.get_mp(), get_mp(), t.get_mp());
    return r;
  }
  BigRat operator+() const {
    return BigRat(*this);
  }
  BigRat operator-() const {
    BigRat r;
    mpq_neg(r.get_mp(), get_mp());
    return r;
  }
  BigRat& operator++() {
    makeCopy();
    mpz_add(get_num_mp(),get_num_mp(),get_den_mp());
    return *this;
  }
  BigRat& operator--() {
    makeCopy();
    mpz_sub(get_num_mp(),get_num_mp(),get_den_mp());
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
    mpq_canonicalize(get_mp());
  } 
  /// Has Exact Division
  static bool hasExactDivision() {
    return true;
  }

  /// return mpz pointer of numerator (const)
  mpz_srcptr get_num_mp() const {
    return mpq_numref(get_mp());
  }
  /// return mpz pointer of numerator
  mpz_ptr get_num_mp() {
    return mpq_numref(get_mp());
  }
  /// return mpz pointer of denominator
  mpz_srcptr get_den_mp() const {
    return mpq_denref(get_mp());
  }
  /// return mpz pointer of denominator
  mpz_ptr get_den_mp() {
    return mpq_denref(get_mp());
  }

  /// get mpq pointer (const)
  mpq_srcptr get_mp() const {
    return rep->get_mp();
  }
  /// get mpq pointer
  mpq_ptr get_mp() {
    return rep->get_mp();
  }
  //@}

  /// \name String Conversion Functions
  //@{
  /// set value from <tt>const char*</tt>
  int set_str(const char* s, int base = 0) {
    makeCopy();
    return mpq_set_str(get_mp(), s, base);
  }
  /// convert to <tt>std::string</tt>
  std::string get_str(int base = 10) const {
    int n = mpz_sizeinbase(mpq_numref(get_mp()), base) + mpz_sizeinbase(mpq_denref(get_mp()), base)+ 3;
    char *buffer = new char[n];
    mpq_get_str(buffer, base, get_mp());
    std::string result(buffer);
    delete [] buffer;
    return result;
  }
  //@}

  /// \name Conversion Functions
  //@{
  /// intValue
  int intValue() const {
    return static_cast<int>(doubleValue());
  }
  /// longValue
  long longValue() const {
    return static_cast<long>(doubleValue());
  }
  /// doubleValue
  double doubleValue() const {
    return mpq_get_d(get_mp());
  }
  /// BigIntValue
  BigInt BigIntValue() const {
    BigInt r;
    mpz_tdiv_q(r.get_mp(), get_num_mp(), get_den_mp());
    return r;
  }
  //@}
}; //BigRat class

inline BigRat operator+(const BigRat& a, const BigRat& b) {
  BigRat r;
  mpq_add(r.get_mp(), a.get_mp(), b.get_mp());
  return r;
}
inline BigRat operator-(const BigRat& a, const BigRat& b) {
  BigRat r;
  mpq_sub(r.get_mp(), a.get_mp(), b.get_mp());
  return r;
}
inline BigRat operator*(const BigRat& a, const BigRat& b) {
  BigRat r;
  mpq_mul(r.get_mp(), a.get_mp(), b.get_mp());
  return r;
}
inline BigRat operator/(const BigRat& a, const BigRat& b) {
  BigRat r;
  mpq_div(r.get_mp(), a.get_mp(), b.get_mp());
  return r;
}
// Chee (3/19/2004):
//   The following definitions of div_exact(x,y) and gcd(x,y)
//   ensures that in Polynomial<NT>
/// divisible(x,y) = "x | y"
inline BigRat div_exact(const BigRat& x, const BigRat& y) {
	BigRat z;
	mpq_div(z.get_mp(), x.get_mp(), y.get_mp());
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
  //	return BigRat(1);  // Remark: we may want replace this by
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
  mpq_div(r.get_mp(), x.get_mp(), y.get_mp());
  return isInteger(r);
}
inline BigRat operator<<(const BigRat& a, unsigned long ul) {
  BigRat r;
  mpq_mul_2exp(r.get_mp(), a.get_mp(), ul);
  return r;
}
inline BigRat operator>>(const BigRat& a, unsigned long ul) {
  BigRat r;
  mpq_div_2exp(r.get_mp(), a.get_mp(), ul);
  return r;
}

inline int cmp(const BigRat& x, const BigRat& y) {
  return mpq_cmp(x.get_mp(), y.get_mp());
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
  //return CORE::operator<<(o, x.get_mp());
  return CORE::io_write(o, x.get_mp());
}
inline std::istream& operator>>(std::istream& i, BigRat& x) {
  x.makeCopy();
  //return CORE::operator>>(i, x.get_mp());
  return CORE::io_read(i, x.get_mp());
}

/// sign
inline int sign(const BigRat& a) {
  return mpq_sgn(a.get_mp());
}
/// abs
inline BigRat abs(const BigRat& a) {
  BigRat r;
  mpq_abs(r.get_mp(), a.get_mp());
  return r;
}
/// neg
inline BigRat neg(const BigRat& a) {
  BigRat r;
  mpq_neg(r.get_mp(), a.get_mp());
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
