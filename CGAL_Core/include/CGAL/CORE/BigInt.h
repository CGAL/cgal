/****************************************************************************
 * Core Library Version 1.7, August 2004
 * Copyright (c) 1995-2004 Exact Computation Project
 * All rights reserved.
 *
 * This file is part of CGAL (www.cgal.org).
 *
 * File: BigInt.h
 * Synopsis:
 *                 a wrapper class for mpz from GMP
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
#ifndef _CORE_BIGINT_H_
#define _CORE_BIGINT_H_

#include <CGAL/CORE/Gmp.h>
#include <CGAL/CORE/RefCount.h>
#include <CGAL/CORE/MemoryPool.h>
#include <string>

namespace CORE {


class BigIntRep : public RCRepImpl<BigIntRep> {
public:
  BigIntRep() {
    mpz_init(mp);
  }
  // Note : should the copy-ctor be alloed at all ? [Sylvain Pion]
  BigIntRep(const BigIntRep& z) : RCRepImpl<BigIntRep>() {
    mpz_init_set(mp, z.mp);
  }
  BigIntRep(signed char c) {
    mpz_init_set_si(mp, c);
  }
  BigIntRep(unsigned char c) {
    mpz_init_set_ui(mp, c);
  }
  BigIntRep(signed int i) {
    mpz_init_set_si(mp, i);
  }
  BigIntRep(unsigned int i) {
    mpz_init_set_ui(mp, i);
  }
  BigIntRep(signed short int s) {
    mpz_init_set_si(mp, s);
  }
  BigIntRep(unsigned short int s) {
    mpz_init_set_ui(mp, s);
  }
  BigIntRep(signed long int l) {
    mpz_init_set_si(mp, l);
  }
  BigIntRep(unsigned long int l) {
    mpz_init_set_ui(mp, l);
  }
  BigIntRep(float f) {
    mpz_init_set_d(mp, f);
  }
  BigIntRep(double d) {
    mpz_init_set_d(mp, d);
  }
  BigIntRep(const char* s, int base=0) {
    mpz_init_set_str(mp, s, base);
  }
  BigIntRep(const std::string& s, int base=0) {
    mpz_init_set_str(mp, s.c_str(), base);
  }
  explicit BigIntRep(mpz_srcptr z) {
    mpz_init_set(mp, z);
  }
  ~BigIntRep() {
    mpz_clear(mp);
  }

  CGAL_CORE_EXPORT CORE_NEW(BigIntRep)
  CGAL_CORE_EXPORT CORE_DELETE(BigIntRep)

  mpz_srcptr get_mp() const {
    return mp;
  }
  mpz_ptr get_mp() {
    return mp;
  }
private:
  mpz_t mp;
};

typedef RCImpl<BigIntRep> RCBigInt;
class CGAL_CORE_EXPORT BigInt : public RCBigInt {
public:
  /// \name Constructors
  //@{
  /// default constructor
  BigInt() : RCBigInt(new BigIntRep()) {}
  /// constructor for <tt>signed char</tt>
  BigInt(signed char x) : RCBigInt(new BigIntRep(x)) {}
  /// constructor for <tt>unsigned char</tt>
  BigInt(unsigned char x) : RCBigInt(new BigIntRep(x)) {}
  /// constructor for <tt>signed short int</tt>
  BigInt(signed short int x) : RCBigInt(new BigIntRep(x)) {}
  /// constructor for <tt>unsigned short int</tt>
  BigInt(unsigned short int x) : RCBigInt(new BigIntRep(x)) {}
  /// constructor for <tt>signed int</tt>
  BigInt(signed int x) : RCBigInt(new BigIntRep(x)) {}
  /// constructor for <tt>unsigned int</tt>
  BigInt(unsigned int x) : RCBigInt(new BigIntRep(x)) {}
  /// constructor for <tt>signed long int</tt>
  BigInt(signed long int x) : RCBigInt(new BigIntRep(x)) {}
  /// constructor for <tt>unsigned long int</tt>
  BigInt(unsigned long int x) : RCBigInt(new BigIntRep(x)) {}
  /// constructor for <tt>float</tt>
  BigInt(float x) : RCBigInt(new BigIntRep(x)) {}
  /// constructor for <tt>double</tt>
  BigInt(double x) : RCBigInt(new BigIntRep(x)) {}
  /// constructor for <tt>const char*</tt> with base
  BigInt(const char* s, int base=0) : RCBigInt(new BigIntRep(s, base)) {}
  /// constructor for <tt>std::string</tt> with base
  BigInt(const std::string& s, int base=0) : RCBigInt(new BigIntRep(s, base)) {}
  /// constructor for <tt>mpz_srcptr</tt>
  explicit BigInt(mpz_srcptr z) : RCBigInt(new BigIntRep(z)) {}
  //@}

  /// \name Copy-Assignment-Destructor
  //@{
  /// copy constructor
  BigInt(const BigInt& rhs) : RCBigInt(rhs) {
    rep->incRef();
  }
  /// assignment operator
  BigInt& operator=(const BigInt& rhs) {
    if (this != &rhs) {
      rep->decRef();
      rep = rhs.rep;
      rep->incRef();
    }
    return *this;
  }
  /// destructor
  ~BigInt() {
    rep->decRef();
  }
  //@}

  /// \name Overloaded operators
  //@{
  BigInt& operator +=(const BigInt& rhs) {
    makeCopy();
    mpz_add(get_mp(), get_mp(), rhs.get_mp());
    return *this;
  }
  BigInt& operator -=(const BigInt& rhs) {
    makeCopy();
    mpz_sub(get_mp(), get_mp(), rhs.get_mp());
    return *this;
  }
  BigInt& operator *=(const BigInt& rhs) {
    makeCopy();
    mpz_mul(get_mp(), get_mp(), rhs.get_mp());
    return *this;
  }
  BigInt& operator /=(const BigInt& rhs) {
    makeCopy();
    mpz_tdiv_q(get_mp(), get_mp(), rhs.get_mp());
    return *this;
  }
  BigInt& operator %=(const BigInt& rhs) {
    makeCopy();
    mpz_tdiv_r(get_mp(), get_mp(), rhs.get_mp());
    return *this;
  }
  BigInt& operator &=(const BigInt& rhs) {
    makeCopy();
    mpz_and(get_mp(), get_mp(), rhs.get_mp());
    return *this;
  }
  BigInt& operator |=(const BigInt& rhs) {
    makeCopy();
    mpz_ior(get_mp(), get_mp(), rhs.get_mp());
    return *this;
  }
  BigInt& operator ^=(const BigInt& rhs) {
    makeCopy();
    mpz_xor(get_mp(), get_mp(), rhs.get_mp());
    return *this;
  }
  BigInt& operator <<=(unsigned long ul) {
    makeCopy();
    mpz_mul_2exp(get_mp(), get_mp(), ul);
    return *this;
  }
  BigInt& operator >>=(unsigned long ul) {
    makeCopy();
    mpz_tdiv_q_2exp(get_mp(), get_mp(), ul);
    return *this;
  }
  //@}

  /// \name unary, increment, decrement operators
  //@{
  BigInt operator+() const {
    return BigInt(*this);
  }
  BigInt operator-() const {
    BigInt r;
    mpz_neg(r.get_mp(), get_mp());
    return r;
  }
  BigInt& operator++() {
    makeCopy();
    mpz_add_ui(get_mp(), get_mp(), 1);
    return *this;
  }
  BigInt& operator--() {
    makeCopy();
    mpz_sub_ui(get_mp(), get_mp(), 1);
    return *this;
  }
  BigInt operator++(int) {
    BigInt r(*this);
    ++(*this);
    return r;
  }
  BigInt operator--(int) {
    BigInt r(*this);
    --(*this);
    return r;
  }
  //@}

  /// \name Helper Functions
  //@{
  /// Has Exact Division
  static bool hasExactDivision() {
    return false;
  }
  /// get mpz pointer (const)
  mpz_srcptr get_mp() const {
    return rep->get_mp();
  }
  /// get mpz pointer
  mpz_ptr get_mp() {
    return rep->get_mp();
  }
  //@}

  /// \name String Conversion Functions
  //@{
  /// set value from <tt>const char*</tt>
  int set_str(const char* s, int base = 0) {
    makeCopy();
    return mpz_set_str(get_mp(), s, base);
  }
  /// convert to <tt>std::string</tt>
  std::string get_str(int base = 10) const {
    int n = mpz_sizeinbase (get_mp(), base) + 2;
    char *buffer = new char[n];
    mpz_get_str(buffer, base, get_mp());
    std::string result(buffer);
    delete [] buffer;
    return result;
  }
  //@}

  /// \name Conversion Functions
  //@{
  /// intValue
  int intValue() const {
    return static_cast<int>(mpz_get_si(get_mp()));
  }
  /// longValue
  long longValue() const {
    return mpz_get_si(get_mp());
  }
  /// ulongValue
  unsigned long ulongValue() const {
    return mpz_get_ui(get_mp());
  }
  /// doubleValue
  double doubleValue() const {
    return mpz_get_d(get_mp());
  }
  //@}
};

inline BigInt operator+(const BigInt& a, const BigInt& b) {
  BigInt r;
  mpz_add(r.get_mp(), a.get_mp(), b.get_mp());
  return r;
}
inline BigInt operator-(const BigInt& a, const BigInt& b) {
  BigInt r;
  mpz_sub(r.get_mp(), a.get_mp(), b.get_mp());
  return r;
}
inline BigInt operator*(const BigInt& a, const BigInt& b) {
  BigInt r;
  mpz_mul(r.get_mp(), a.get_mp(), b.get_mp());
  return r;
}
inline BigInt operator/(const BigInt& a, const BigInt& b) {
  BigInt r;
  mpz_tdiv_q(r.get_mp(), a.get_mp(), b.get_mp());
  return r;
}
inline BigInt operator%(const BigInt& a, const BigInt& b) {
  BigInt r;
  mpz_tdiv_r(r.get_mp(), a.get_mp(), b.get_mp());
  return r;
}
inline BigInt operator&(const BigInt& a, const BigInt& b) {
  BigInt r;
  mpz_and(r.get_mp(), a.get_mp(), b.get_mp());
  return r;
}
inline BigInt operator|(const BigInt& a, const BigInt& b) {
  BigInt r;
  mpz_ior(r.get_mp(), a.get_mp(), b.get_mp());
  return r;
}
inline BigInt operator^(const BigInt& a, const BigInt& b) {
  BigInt r;
  mpz_xor(r.get_mp(), a.get_mp(), b.get_mp());
  return r;
}
inline BigInt operator<<(const BigInt& a, unsigned long ul) {
  BigInt r;
  mpz_mul_2exp(r.get_mp(), a.get_mp(), ul);
  return r;
}
inline BigInt operator>>(const BigInt& a, unsigned long ul) {
  BigInt r;
  mpz_tdiv_q_2exp(r.get_mp(), a.get_mp(), ul);
  return r;
}

inline int cmp(const BigInt& x, const BigInt& y) {
  return mpz_cmp(x.get_mp(), y.get_mp());
}
inline bool operator==(const BigInt& a, const BigInt& b) {
  return cmp(a, b) == 0;
}
inline bool operator!=(const BigInt& a, const BigInt& b) {
  return cmp(a, b) != 0;
}
inline bool operator>=(const BigInt& a, const BigInt& b) {
  return cmp(a, b) >= 0;
}
inline bool operator>(const BigInt& a, const BigInt& b) {
  return cmp(a, b) > 0;
}
inline bool operator<=(const BigInt& a, const BigInt& b) {
  return cmp(a, b) <= 0;
}
inline bool operator<(const BigInt& a, const BigInt& b) {
  return cmp(a, b) < 0;
}

inline std::ostream& operator<<(std::ostream& o, const BigInt& x) {
  //return CORE::operator<<(o, x.get_mp());
  return CORE::io_write(o, x.get_mp());
}
inline std::istream& operator>>(std::istream& i, BigInt& x) {
  x.makeCopy();
  //return CORE::operator>>(i, x.get_mp());
  return CORE::io_read(i, x.get_mp());
}

/// sign
inline int sign(const BigInt& a) {
  return mpz_sgn(a.get_mp());
}
/// abs
inline BigInt abs(const BigInt& a) {
  BigInt r;
  mpz_abs(r.get_mp(), a.get_mp());
  return r;
}
/// neg
inline BigInt neg(const BigInt& a) {
  BigInt r;
  mpz_neg(r.get_mp(), a.get_mp());
  return r;
}
/// negate
inline void negate(BigInt& a) {
  a.makeCopy();
  mpz_neg(a.get_mp(), a.get_mp());
}
/// cmpabs
inline int cmpabs(const BigInt& a, const BigInt& b) {
  return mpz_cmpabs(a.get_mp(), b.get_mp());
}

/// \name Conversion Functions
//@{
/// longValue
inline long longValue(const BigInt& a) {
  return a.longValue();
}
/// ulongValue
inline unsigned long ulongValue(const BigInt& a) {
  return a.ulongValue();
}
/// doubleValue
inline double doubleValue(const BigInt& a) {
  return a.doubleValue();
}
//@}

/// \name File I/O Functions
//@{
/// read from file
void readFromFile(BigInt& z, std::istream& in, long maxLength = 0);
/// write to file
void writeToFile(const BigInt& z, std::ostream& in, int base=10, int charsPerLine=80);
//@}

/// \name Misc Functions
//@{
/// isEven
inline bool isEven(const BigInt& z) {
  return mpz_even_p(z.get_mp());
}
/// isOdd
inline bool isOdd(const BigInt& z) {
  return mpz_odd_p(z.get_mp());
}

/// get exponent of power 2
inline unsigned long getBinExpo(const BigInt& z) {
  return mpz_scan1(z.get_mp(), 0);
}
/// get exponent of power k
inline void getKaryExpo(const BigInt& z, BigInt& m, int& e, unsigned long k) {
  mpz_t f;
  mpz_init_set_ui(f, k);
  m.makeCopy();
  e = mpz_remove(m.get_mp(), z.get_mp(), f);
  mpz_clear(f);
}

/// divisible(x,y) = "x | y"
inline bool isDivisible(const BigInt& x, const BigInt& y) {
  return mpz_divisible_p(x.get_mp(), y.get_mp()) != 0;
}
inline bool isDivisible(int x, int y) {
  return x % y == 0;
}
inline bool isDivisible(long x, long y) {
  return x % y == 0;
}
/// exact div
inline void divexact(BigInt& z, const BigInt& x, const BigInt& y) {
  z.makeCopy();
  mpz_divexact(z.get_mp(), x.get_mp(), y.get_mp());
}
// Chee (1/12/2004)   The definition of div_exact(x,y) next
//   ensure that in Polynomials<NT> works with both NT=BigInt and NT=int:
inline BigInt div_exact(const BigInt& x, const BigInt& y) {
  BigInt z;             // precodition: isDivisible(x,y)
  divexact(z, x, y); // z is set to x/y;
  return z;
}
inline int div_exact(int x, int y) {
  return x/y;  // precondition: isDivisible(x,y)
}
inline long div_exact(long x, long y) {
  return x/y;  // precondition: isDivisible(x,y)
}


/// gcd
inline BigInt gcd(const BigInt& a, const BigInt& b) {
  BigInt r;
  mpz_gcd(r.get_mp(), a.get_mp(), b.get_mp());
  return r;
}
/// div_rem
inline void div_rem(BigInt& q, BigInt& r, const BigInt& a, const BigInt& b) {
  q.makeCopy();
  r.makeCopy();
  mpz_tdiv_qr(q.get_mp(), r.get_mp(), a.get_mp(), b.get_mp());
}
/// power
inline void power(BigInt& c, const BigInt& a, unsigned long ul) {
  c.makeCopy();
  mpz_pow_ui(c.get_mp(), a.get_mp(), ul);
}

// pow
inline BigInt pow(const BigInt& a, unsigned long ui) {
  BigInt r;
  power(r, a, ui);
  return r;
}

// bit length
inline int bitLength(const BigInt& a) {
  return mpz_sizeinbase(a.get_mp(), 2);
}
/// floorLg -- floor of log_2(a)
/** Convention: a=0, floorLg(a) returns -1.
 *  This makes sense for integer a.
 */
inline long floorLg(const BigInt& a) {
  return (sign(a) == 0) ? (-1) : (bitLength(a)-1);
}
/// ceilLg -- ceiling of log_2(a) where a=BigInt, int or long
/** Convention: a=0, ceilLg(a) returns -1.
 *  This makes sense for integer a.
 */
inline long ceilLg(const BigInt& a) {
  if (sign(a) == 0)
    return -1;
  unsigned long len = bitLength(a);
  return (mpz_scan1(a.get_mp(), 0) == len-1) ? (len-1) : len;
}
inline long ceilLg(long a) { // need this for Polynomial<long>
  return ceilLg(BigInt(a));
}
inline long ceilLg(int a) { // need this for Polynomial<int>
  return ceilLg(BigInt(a));
}

//@}



} //namespace CORE

#endif // _CORE_BIGINT_H_
