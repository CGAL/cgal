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
 * File: BigFloat.h
 * Synopsis: 
 * 		An implementation of BigFloat numbers with error bounds.
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

#ifndef _CORE_BIGFLOAT_H_
#define _CORE_BIGFLOAT_H_

#include <CGAL/CORE/BigFloatRep.h>

namespace CORE { 

class Expr;

/// \class BigFloat BigFloat.h
/// \brief BigFloat is a class of Float-Point number with error bounds.
typedef RCImpl<BigFloatRep> RCBigFloat;

class CGAL_CORE_EXPORT BigFloat : public RCBigFloat {
public:
  /// \name Constructors and Destructor
  //@{
  /// default constructor
  BigFloat() : RCBigFloat(new BigFloatRep()) {}
  /// constructor for <tt>short</tt>
  BigFloat(short i) : RCBigFloat(new BigFloatRep(i)) {}
  /// constructor for <tt>float</tt>
  BigFloat(float i) : RCBigFloat(new BigFloatRep(i)) {}
  /// constructor for <tt>int</tt>
  BigFloat(int i) : RCBigFloat(new BigFloatRep(i)) {}
  /// constructor for <tt>long</tt>
  BigFloat(long l) : RCBigFloat(new BigFloatRep(l)) {}
  /// constructor for <tt>double</tt>
  BigFloat(double d) : RCBigFloat(new BigFloatRep(d)) {}
  /// constructor for <tt>const char* </tt>(default base = 10)
  BigFloat(const char* s) : RCBigFloat(new BigFloatRep(s)) {}
  /// constructor for <tt>std::string</tt>(default base = 10)
  BigFloat(const std::string& s) : RCBigFloat(new BigFloatRep(s)) {}

  /// constructor for <tt>int</tt> and <tt>long</tt>
  //     This is a hack because in Sturm, we need to approximate any
  //     coefficient type NT to a BigFloat, and it would complain if we
  //     do not have this method explicitly:
  BigFloat(int& i, const extLong& /*r*/, const extLong& /*a*/)
      : RCBigFloat(new BigFloatRep(i)) {}
  BigFloat(long& x, const extLong& /*r*/, const extLong& /*a*/)
      : RCBigFloat(new BigFloatRep(x)) {}
  /// constructor from <tt>BigInt</tt>, error and exponent values
  BigFloat(const BigInt& I, unsigned long er, long ex)
      : RCBigFloat(new BigFloatRep(I, er, ex)) {}
  /// constructor from <tt>BigInt</tt>, exponent values
  BigFloat(const BigInt& I, long ex)
      : RCBigFloat(new BigFloatRep(I, ex)) {}
  BigFloat(const BigInt& I)
      : RCBigFloat(new BigFloatRep(I)) {}
  /// constructor for <tt>BigRat</tt>
  BigFloat(const BigRat& R, const extLong& r = defRelPrec,
           const extLong& a = defAbsPrec)
      : RCBigFloat(new BigFloatRep()) {
    rep->approx(R, r, a);
  }

  // REMARK: it is somewhat against our principles to have BigFloat
  // know about Expr, but BigFloat has a special role in our system!
  // ===============================
  /// constructor for <tt>Expr</tt>
  explicit BigFloat(const Expr& E, const extLong& r = defRelPrec,
           const extLong& a = defAbsPrec);

  //Dummy
  explicit BigFloat(const BigFloat& E, const extLong& ,
           const extLong&): RCBigFloat(E) {
    rep->incRef();
  }

  /// constructor for <tt>BigFloatRep</tt>
  explicit BigFloat(BigFloatRep* r) : RCBigFloat(new BigFloatRep()) {
    rep = r;
  }
  //@}

  /// \name Copy-Assignment-Destructor
  //@{
  /// copy constructor
  BigFloat(const BigFloat& rhs) : RCBigFloat(rhs) {
    rep->incRef();
  }
  /// assignment operator
  BigFloat& operator=(const BigFloat& rhs) {
    if (this != &rhs) {
      rep->decRef();
      rep = rhs.rep;
      rep->incRef();
    }
    return *this;
  }
  /// destructor
  ~BigFloat() {
    rep->decRef();
  }
  //@}

  /// \name Compound Assignment Operators
  //@{
  /// operator+=
  BigFloat& operator+= (const BigFloat& x) {
    BigFloat z;
    z.rep->add(*rep, *x.rep);
    *this = z;
    return *this;
  }
  /// operator-=
  BigFloat& operator-= (const BigFloat& x) {
    BigFloat z;
    z.rep->sub(*rep, *x.rep);
    *this = z;
    return *this;
  }
  /// operator*=
  BigFloat& operator*= (const BigFloat& x) {
    BigFloat z;
    z.rep->mul(*rep, *x.rep);
    *this = z;
    return *this;
  }
  /// operator/=
  BigFloat& operator/= (const BigFloat& x) {
    BigFloat z;
    z.rep->div(*rep, *x.rep, defBFdivRelPrec);
    *this = z;
    return *this;
  }
  //@}

  /// \name Unary Minus Operator
  //@{
  /// unary plus
  BigFloat operator+() const {
    return BigFloat(*this);
  }
  /// unary minus
  BigFloat operator-() const {
    return BigFloat(-rep->m, rep->err, rep->exp);
  }
  //@}

  /// \name String Conversion Functions
  //@{
  /// set value from <tt>const char*</tt> (base = 10)
  void fromString(const char* s, const extLong& p=defBigFloatInputDigits) {
    rep->fromString(s, p);
  }
  /// convert to <tt>std::string</tt> (base = 10)
  std::string toString(long prec=defBigFloatOutputDigits, bool sci=false) const {
    return rep->toString(prec, sci);
  }
  std::string str() const {
    return toString();
  }
  //@}

  /// \name Conversion Functions
  //@{
  /// return int value
  int intValue() const {
    return static_cast<int>(rep->toLong());
  }
  /// return long value
  long longValue() const {
    long l = rep->toLong();
    if ((l == LONG_MAX) || (l == LONG_MIN))
      return l; // return the overflown value.
    if ((sign() < 0) && (cmp(BigFloat(l)) != 0)) {
      // a negative value not exactly rounded.
      l--; // rounded to floor.
    }
    return l;
  }
  /// return float value
  float floatValue() const {
    return static_cast<float>(rep->toDouble());
  }
  /// return double value
  double doubleValue() const {
    return rep->toDouble();
  }
  /// return BigInt value
  BigInt BigIntValue() const {
    return rep->toBigInt();
  }
  /// return BigRat value
  BigRat BigRatValue() const {
    return rep->BigRatize();
  }
  //@}

  /// \name Helper Functions
  //@{
  /// Has Exact Division
  static bool hasExactDivision() {
    return false;
  }
  
  //CONSTANTS
  /// return BigFloat(0)
  static const BigFloat& getZero();
  /// return BigFloat(1)
  static const BigFloat& getOne();

  /// sign function
  /** \note This is only the sign of the mantissa, it can be taken to be
      the sign of the BigFloat only if !(isZeroIn()). */
  int sign() const {
    assert((err() == 0 && m() == 0) || !(isZeroIn()));
    return rep->signM();
  }
  /// check whether contains zero
  /** \return true if contains zero, otherwise false */
  bool isZeroIn() const {
    return rep->isZeroIn();
  }
  /// absolute value function
  BigFloat abs() const {
    return (sign()>0) ? +(*this) : -(*this);
  }
  ///  comparison function
  int cmp(const BigFloat& x) const {
    return rep->compareMExp(*x.rep);
  }

  /// get mantissa
  const BigInt& m() const {
    return rep->m;
  }
  /// get error bits
  unsigned long err() const {
    return rep->err;
  }
  /// get exponent
  long exp() const {
    return rep->exp;
  }

  /// check whether err == 0
  /** \return true if err == 0, otherwise false */
  bool isExact() const {
    return rep->err == 0;
  }
  /// set err to 0
  /** \return an exact BigFloat, see Tutorial for why this is useful! */
  BigFloat& makeExact() {
    makeCopy();
    rep->err =0;
    return *this;
  }
  /// set err to 0, but first add err to the mantissa (m)
  /** \return the ceiling exact BigFloat, variant of makeExact */
  BigFloat& makeCeilExact() {
    makeCopy();
    rep->m += rep->err;
    rep->err =0;
    return *this;
  }
  /// set err to 0, but subtract err from the mantissa (m)
  /** \return the floor exact BigFloat, variant of makeExact */
  BigFloat& makeFloorExact() {
    makeCopy();
    rep->m -= rep->err;
    rep->err =0;
    return *this;
  }
  /// set err to 1
  /** \return an inexact BigFloat, see Tutorial for why this is useful! */
  BigFloat& makeInexact() {
    makeCopy();
    rep->err =1;
    return *this;
  }

  /// return lower bound of Most Significant Bit
  extLong lMSB() const {
    return rep->lMSB();
  }
  /// return upper bound of Most Significant Bit
  extLong uMSB() const {
    return rep->uMSB();
  }
  /// return Most Significant Bit
  extLong MSB() const {
    return rep->MSB();
  }

  /// floor of Lg(err)
  extLong flrLgErr() const {
    return rep->flrLgErr();
  }
  /// ceil of Lg(err)
  extLong clLgErr() const {
    return rep->clLgErr();
  }

  /// division with relative precsion <tt>r</tt>
  BigFloat div(const BigFloat& x, const extLong& r) const {
    BigFloat y;
    y.rep->div(*rep, *x.rep, r);
    return y;
  }
  /// exact division by 2
  BigFloat div2() const {
    BigFloat y;
    y.rep->div2(*rep);
    return y;
  }

  /// squareroot
  BigFloat sqrt(const extLong& a) const {
    BigFloat x;
    x.rep->sqrt(*rep, a);
    return x;
  }
  /// squareroot with initial approximation <tt>init</tt>
  BigFloat sqrt(const extLong& a, const BigFloat& init) const {
    BigFloat x;
    x.rep->sqrt(*rep, a, init);
    return x;
  }
  //@}

  /// \name Utility Functions
  //@{
  /// approximate BigInt number
  void approx(const BigInt& I, const extLong& r, const extLong& a) {
    makeCopy();
    rep->trunc(I, r, a);
  }
  /// approximate BigFloat number
  void approx(const BigFloat& B, const extLong& r, const extLong& a) {
    makeCopy();
    rep->approx(*B.rep, r, a);
  }
  /// approximate BigRat number
  void approx(const BigRat& R, const extLong& r, const extLong& a) {
    makeCopy();
    rep->approx(R, r, a);
  }
  /// dump internal data
  void dump() const {
    rep->dump();
  }
  //@}

  /// returns a BigFloat of value \f$ 2^e \f$
  static BigFloat exp2(int e) {
    return BigFloat(BigFloatRep::exp2(e));
  }

}; // class BigFloat

//@} // For compatibility with BigInt

/// \name File I/O Functions
//@{
/// read from file
void readFromFile(BigFloat& bf, std::istream& in, long maxLength = 0);
/// write to file
void writeToFile(const BigFloat& bf, std::ostream& in,
		int base=10, int charsPerLine=80);

/// IO stream operator<<
inline std::ostream& operator<< (std::ostream& o, const BigFloat& x) {
  x.getRep().operator<<(o);
  return o;
}
/// IO stream operator>>
inline std::istream& operator>> (std::istream& i, BigFloat& x) {
  x.makeCopy();
  x.getRep().operator>>(i);
  return i;
}
//@}

/// operator+
inline BigFloat operator+ (const BigFloat& x, const BigFloat& y) {
  BigFloat z;
  z.getRep().add(x.getRep(), y.getRep());
  return z;
}
/// operator-
inline BigFloat operator- (const BigFloat& x, const BigFloat& y) {
  BigFloat z;
  z.getRep().sub(x.getRep(), y.getRep());
  return z;
}
/// operator*
inline BigFloat operator* (const BigFloat& x, const BigFloat& y) {
  BigFloat z;
  z.getRep().mul(x.getRep(), y.getRep());
  return z;
}
/// operator/
inline BigFloat operator/ (const BigFloat& x, const BigFloat& y) {
  BigFloat z;
  z.getRep().div(x.getRep(),y.getRep(),defBFdivRelPrec);
  return z;
}

/// operator==
inline bool operator== (const BigFloat& x, const BigFloat& y) {
  return x.cmp(y) == 0;
}
/// operator!=
inline bool operator!= (const BigFloat& x, const BigFloat& y) {
  return x.cmp(y) != 0;
}
/// operator>=
inline bool operator>= (const BigFloat& x, const BigFloat& y) {
  return x.cmp(y) >= 0;
}
/// operator>
inline bool operator> (const BigFloat& x, const BigFloat& y) {
  return x.cmp(y) > 0;
}
/// operator<=
inline bool operator<= (const BigFloat& x, const BigFloat& y) {
  return x.cmp(y) <= 0;
}
/// operator<
inline bool operator< (const BigFloat& x, const BigFloat& y) {
  return x.cmp(y) < 0;
}

/// sign
inline int sign(const BigFloat& x) {
  return x.sign();
}
/// div2
inline BigFloat div2(const BigFloat& x){
  return x.div2();
}
/// abs
inline BigFloat abs(const BigFloat& x) {
  return x.abs();
}
/// cmp
inline int cmp(const BigFloat& x, const BigFloat& y) {
  return x.cmp(y);
}
/// pow
CGAL_CORE_EXPORT BigFloat pow(const BigFloat&, unsigned long);
/// power
inline BigFloat power(const BigFloat& x, unsigned long p) {
  return pow(x, p);
}
/// root(x,k,prec,xx) returns the k-th root of x to precision prec.
///   The argument x is an initial approximation.
BigFloat root(const BigFloat&, unsigned long k, const extLong&, const BigFloat&);
inline BigFloat root(const BigFloat& x, unsigned long k) {
  return root(x, k, defBFsqrtAbsPrec, x);
}

/// sqrt to defAbsPrec:
inline BigFloat sqrt(const BigFloat& x) {
  return x.sqrt(defBFsqrtAbsPrec);
}

/// convert an BigFloat Interval to a BigFloat with error bits
inline BigFloat centerize(const BigFloat& a, const BigFloat& b) {
  BigFloat z;
  z.getRep().centerize(a.getRep(), b.getRep());
  return z;
}

/// minStar(m,n) returns the min-star of m and n
inline long minStar(long m, long n) {
  if (m*n <= 0) return 0;
  if (m>0) 
    return core_min(m, n);
  else 
    return core_max(m, n);
}
/// \name Functions for Compatibility with BigInt (needed by Poly, Curves)
//@{
/// isDivisible(a,b) = "is a divisible by b"
/** 	Assuming that a and  b are in coanonized forms.
	Defined to be true if mantissa(b) | mantissa(a) && 
	exp(b) = min*(exp(b), exp(a)).
 *      This concepts assume a and b are exact BigFloats.
 */
inline bool isDivisible(const BigFloat& a, const BigFloat& b) {
  // assert: a and b are exact BigFloats.
  if (sign(a.m()) == 0) return true;
  if (sign(b.m()) == 0) return false;
  unsigned long bin_a = getBinExpo(a.m());
  unsigned long bin_b = getBinExpo(b.m());
  
  BigInt m_a = a.m() >> bin_a;
  BigInt m_b = b.m() >> bin_b;
  long e_a = bin_a + BigFloatRep::bits(a.exp());
  long e_b = bin_b + BigFloatRep::bits(b.exp());
  long dx = minStar(e_a, e_b);

  return isDivisible(m_a, m_b) && (dx == e_b); 
}

inline bool isDivisible(double x, double y) {
  //Are these exact?
  return isDivisible(BigFloat(x), BigFloat(y)); 
}

/// div_exact(x,y) returns the BigFloat quotient of x divided by y
/**	This is defined only if isDivisible(x,y).
 */
// Chee (8/1/2004)   The definition of div_exact(x,y) 
//   ensure that Polynomials<NT> works with NT=BigFloat and NT=double.
//Bug: We should first normalize the mantissas of the Bigfloats and
//then do the BigInt division. For e.g. 1 can be written as 2^{14}*2{-14}.
//Now if we divide 2 by one using this representation of one and without
// normalizing it then we get zero.
inline BigFloat div_exact(const BigFloat& x, const BigFloat& y) {
  BigInt z;
  assert (isDivisible(x,y));
  unsigned long bin_x = getBinExpo(x.m());
  unsigned long bin_y = getBinExpo(y.m());

  BigInt m_x = x.m() >> bin_x;
  BigInt m_y = y.m() >> bin_y;
  long e_x = bin_x + BigFloatRep::bits(x.exp());
  long e_y = bin_y + BigFloatRep::bits(y.exp());
  //Since y divides x, e_y = minstar(e_x, e_y)
  z = div_exact(m_x, m_y);

  //  mpz_divexact(z.get_mp(), x.m().get_mp(), y.m().get_mp()); THIS WAS THE BUG
  // assert: x.exp() - y.exp() does not under- or over-flow.
  return BigFloat(z, e_x - e_y);  
}

inline BigFloat div_exact(double x, double y) {
  return div_exact(BigFloat(x), BigFloat(y));
}
// Remark: there is another notion of "exact division" for BigFloats,
// 	and that is to make the division return an "exact" BigFloat
// 	i.e., err()=0.  

/// gcd(a,b) =  BigFloat(gcd(a.mantissa,b.matissa), min(a.exp(), b.exp()) )
inline BigFloat gcd(const BigFloat& a, const BigFloat& b) {
  if (sign(a.m()) == 0) return core_abs(b);
  if (sign(b.m()) == 0) return core_abs(a);

  BigInt r;
  long dx;
  unsigned long bin_a = getBinExpo(a.m());
  unsigned long bin_b = getBinExpo(b.m());

/* THE FOLLOWING IS ALTERNATIVE CODE, for GCD using base B=2^{14}:
 *std::cout << "bin_a=" << bin_a << ",bin_b=" << bin_b << std::endl;
  std::cout << "a.exp()=" << a.exp() << ",b.exp()=" << b.exp() << std::endl;
  long chunk_a = BigFloatRep::chunkFloor(bin_a);
  long chunk_b = BigFloatRep::chunkFloor(bin_b);
  BigInt m_a = BigFloatRep::chunkShift(a.m(), chunk_a);
  BigInt m_b = BigFloatRep::chunkShift(b.m(), chunk_b);

  r = gcd(m_a, m_b);
  dx = minStar(chunk_a + a.exp(), chunk_b + b.exp());
*/
  BigInt m_a = a.m() >> bin_a;
  BigInt m_b = b.m() >> bin_b;
  r = gcd(m_a, m_b);
  dx = minStar(bin_a + BigFloatRep::bits(a.exp()),
		  bin_b + BigFloatRep::bits(b.exp()));

  long chunks = BigFloatRep::chunkFloor(dx);
  r <<= (dx - BigFloatRep::bits(chunks));
  dx = chunks;

  return BigFloat(r, 0, dx);
}

// Not needed for now:
/// div_rem
// inline void div_rem(BigFloat& q, BigFloat& r,
// 	const BigFloat& a, const BigFloat& b) {
  //q.makeCopy();
  //r.makeCopy();
  //mpz_tdiv_qr(q.get_mp(), r.get_mp(), a.get_mp(), b.get_mp());
//}//


// constructor BigRat from BigFloat
inline BigRat::BigRat(const BigFloat& f) : RCBigRat(new BigRatRep()){
  *this = f.BigRatValue();
}
} //namespace CORE
#endif // _CORE_BIGFLOAT_H_
