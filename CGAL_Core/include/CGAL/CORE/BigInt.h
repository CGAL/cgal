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

#include <CGAL/boost_mp.h>
#include <CGAL/CORE/Gmp.h>
#include <CGAL/CORE/RefCount.h>
#include <CGAL/CORE/MemoryPool.h>
#include <string>

namespace CORE {

  typedef  boost::multiprecision::cpp_int Z;



  class BigIntRep : public RCRepImpl<BigIntRep > {
public:
  BigIntRep()
    : mp()
    {}


  // Note : should the copy-ctor be allowed at all ? [Sylvain Pion]
    BigIntRep(const BigIntRep& z) : RCRepImpl<BigIntRep >(), mp(z.mp)
    {}

  BigIntRep(signed char c)
    : mp(c)
    {}

  BigIntRep(unsigned char c)
    : mp(c)
    {}

  BigIntRep(signed int i)
    : mp(i)
    {}

  BigIntRep(unsigned int i)
    : mp(i)
    {}


  BigIntRep(signed short int s)
    : mp(s)
    {}
  BigIntRep(unsigned short int s)
    : mp(s)
    {}

  BigIntRep(signed long int l)
    : mp(l)
    {}

  BigIntRep(unsigned long int l)
    : mp(l)
    {}

  BigIntRep(float f)
    : mp(f)
    {}

  BigIntRep(double d)
    : mp(d)
    {}

  BigIntRep(const char* s, int base=0)
    : mp(s)
    {}

  BigIntRep(const std::string& s, int base=0)
    : mp(s)
    {}


  explicit BigIntRep(const Z& z)
      : mp(z)
  {}
  /*
  ~BigIntRep() {
    mpz_clear(mp);
  }
    */

    //CGAL_CORE_EXPORT CORE_NEW(BigIntRep)
    //CGAL_CORE_EXPORT CORE_DELETE(BigIntRep)

  const Z& get_mp() const {
    return mp;
  }
  Z&  get_mp() {
    return mp;
  }
private:
  Z mp;
};

  //typedef RCImpl<BigIntRep> RCBigInt;

class CGAL_CORE_EXPORT BigInt : public RCImpl<BigIntRep> {
public:

    typedef RCImpl<BigIntRep> RCBigInt;

  /// \name Constructors
  //@{
  /// default constructor
  BigInt() : RCBigInt(new BigIntRep()) {}
  BigInt(const Z& z) : RCBigInt(new BigIntRep(z)) {}
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
  // explicit BigInt(mpz_srcptr z) : RCBigInt(new BigIntRep(z)) {}
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
    get_mp() += rhs.get_mp();
    return *this;
  }

  BigInt& operator -=(const BigInt& rhs) {
    makeCopy();
    get_mp() -= rhs.get_mp();
    return *this;
  }
  BigInt& operator *=(const BigInt& rhs) {
    makeCopy();
    get_mp() *= rhs.get_mp();
    return *this;
  }
  BigInt& operator /=(const BigInt& rhs) {
    makeCopy();
    get_mp() /= rhs.get_mp();
    return *this;
  }
  BigInt& operator %=(const BigInt& rhs) {
    makeCopy();
    get_mp() %= rhs.get_mp();
    return *this;
  }
  BigInt& operator &=(const BigInt& rhs) {
    makeCopy();
    get_mp() &= rhs.get_mp();
    return *this;
  }
  BigInt& operator |=(const BigInt& rhs) {
    makeCopy();
    get_mp() |= rhs.get_mp();
    return *this;
  }
  BigInt& operator ^=(const BigInt& rhs) {
    makeCopy();
    get_mp() ^= rhs.get_mp();
    return *this;
  }
  BigInt& operator <<=(unsigned long ul) {
    makeCopy();
    get_mp() <<= ul;
    return *this;
  }
  BigInt& operator >>=(unsigned long ul) {
    makeCopy();
    get_mp() >>= ul;
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
    r.get_mp() = -get_mp();
    return r;
  }
  BigInt& operator++() {
    makeCopy();
    ++get_mp();
    return *this;
  }
  BigInt& operator--() {
    makeCopy();
    --get_mp();
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
  const Z& get_mp() const {
    return rep->get_mp();
  }
  /// get mpz pointer
  Z& get_mp() {
    return rep->get_mp();
  }
  //@}

  /// \name String Conversion Functions
  //@{
  /// set value from <tt>const char*</tt>
  int set_str(const char* s, int base = 0) {
    makeCopy();
    get_mp() = Z(s);
    return 0;  // should be -1 if not correct in the base (we ignore)
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
  /// ulongValue
  unsigned long ulongValue() const {
    return get_mp().convert_to<unsigned long>();
  }
  /// doubleValue
  double doubleValue() const {
    return get_mp().convert_to<double>();
  }
  //@}
};

inline BigInt operator+(const BigInt& a, const BigInt& b) {
  BigInt r;
  r.get_mp() = a.get_mp() + b.get_mp();
  return r;
}

inline BigInt operator-(const BigInt& a, const BigInt& b) {
  BigInt r;
  r.get_mp() = a.get_mp() - b.get_mp();
  return r;
}

inline BigInt operator*(const BigInt& a, const BigInt& b) {
  BigInt r;
  r.get_mp() = a.get_mp() * b.get_mp();
  return r;
}

inline BigInt operator/(const BigInt& a, const BigInt& b) {
  BigInt r;
  r.get_mp() = a.get_mp() / b.get_mp();
  return r;
}

inline BigInt operator%(const BigInt& a, const BigInt& b) {
  BigInt r;
  r.get_mp() = a.get_mp() % b.get_mp();
  return r;
}

inline BigInt operator&(const BigInt& a, const BigInt& b) {
  BigInt r;
  r.get_mp() = a.get_mp() & b.get_mp();
  return r;
}

inline BigInt operator|(const BigInt& a, const BigInt& b) {
  BigInt r;
  r.get_mp() = a.get_mp()| b.get_mp();
  return r;
}

inline BigInt operator^(const BigInt& a, const BigInt& b) {
  BigInt r;
  r.get_mp() = a.get_mp() ^ b.get_mp();
  return r;
}

inline BigInt operator<<(const BigInt& a, unsigned long ul) {
  BigInt r;
  r.get_mp() = a.get_mp() << ul;
  return r;
}

inline BigInt operator>>(const BigInt& a, unsigned long ul) {
  BigInt r;
  r.get_mp() = a.get_mp() >> ul;
  return r;
}


inline int cmp(const BigInt& x, const BigInt& y) {
  return x.get_mp().compare(y.get_mp());
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
  return o <<x.get_mp();
}

inline std::istream& operator>>(std::istream& i, BigInt& x) {
  x.makeCopy();
  return i >> x.get_mp();
}

/// sign
inline int sign(const BigInt& a) {
  return sign(a.get_mp());
}

/// abs
inline BigInt abs(const BigInt& a) {
  BigInt r;
  r.get_mp() = abs(a.get_mp());
  return r;
}

/// neg
inline BigInt neg(const BigInt& a) {
  BigInt r;
  r.get_mp() = - a.get_mp();
  return r;
}

/// negate
inline void negate(BigInt& a) {
  a.makeCopy();
  a.get_mp() = - a.get_mp();
}

/// cmpabs
inline int cmpabs(const BigInt& a, const BigInt& b) {
  assert(false);
  // return mpz_cmpabs(a.get_mp(), b.get_mp());   AF: todo
  return 0;
}

/// \name Conversion Functions
//@{
/// longValue
inline long longValue(const BigInt& a) {
  return a.longValue();
}

/// ulongValue
inline unsigned long ulongValue(const BigInt& a) {
    assert(a >= BigInt(0));
  return a.ulongValue();
}

/// doubleValue
inline double doubleValue(const BigInt& a) {
  return a.doubleValue();
}
//@}

/*

/// \name File I/O Functions
//@{
/// read from file

void readFromFile(BigInt& z, std::istream& in, long maxLength = 0);
/// write to file
void writeToFile(const BigInt& z, std::ostream& in, int base=10, int charsPerLine=80);
//@}
*/


/// \name Misc Functions
//@{
/// isEven
inline bool isEven(const BigInt& z) {
  return bit_test(z.get_mp(),0) == 0;
}
/// isOdd
inline bool isOdd(const BigInt& z) {
  return bit_test(z.get_mp(),0) == 1;
}

/// get exponent of power 2
inline unsigned long getBinExpo(const BigInt& z) {
    if (z.get_mp().is_zero()) {
        return (std::numeric_limits<unsigned long>::max)();
    }
    return lsb(abs(z.get_mp()));
}

/// get exponent of power k
inline void getKaryExpo(const BigInt& z, BigInt& m, int& e, unsigned long uk) {
    BigInt k(uk), q, r;
    e = 0;
    m = z;
    for(;;) {
        divide_qr(m.get_mp(), k.get_mp(), q.get_mp(), r.get_mp());
        if (!r.get_mp().is_zero()) break;
        m.get_mp() = q.get_mp();
        ++e;
    }
}

/// divisible(x,y) = "x | y"
inline bool isDivisible(const BigInt& x, const BigInt& y) {
    BigInt z, r;
    divide_qr(x.get_mp(), y.get_mp(), z.get_mp(), r.get_mp());
    return r.get_mp().is_zero();
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
  BigInt r;
  divide_qr(x.get_mp(), y.get_mp(), z.get_mp(), r.get_mp() );  // was void mpz_divexact (mpz_t q, const mpz_t n, const mpz_t d)   Is this faster?
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
  r.get_mp() = gcd(a.get_mp(), b.get_mp());
  return r;
}

/// div_rem
inline void div_rem(BigInt& q, BigInt& r, const BigInt& a, const BigInt& b) {
  q.makeCopy();
  r.makeCopy();
  divide_qr(a.get_mp(), b.get_mp(), q.get_mp(), r.get_mp());
}

/// power
inline void power(BigInt& c, const BigInt& a, unsigned long ul) {
  c.makeCopy();
  c.get_mp() = pow(a.get_mp(), ul);
}

// pow
inline BigInt pow(const BigInt& a, unsigned long ui) {
  BigInt r;
  power(r, a, ui);
  return r;
}

// bit length
inline int bitLength(const BigInt& a) {
    if (a.get_mp().is_zero()) {
        return 0;
    }
  return msb(abs(a.get_mp()))+1;    /// AF todo     was    mpz_sizeinbase(a.get_mp(), 2);
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

  return (lsb(abs(a.get_mp())) == len - 1) ? (len - 1) : len;
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
