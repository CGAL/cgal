/****************************************************************************
 * Core Library Version 1.7, August 2004
 * Copyright (c) 1995-2004 Exact Computation Project
 * All rights reserved.
 *
 * This file is part of CGAL (www.cgal.org);
 *
 * File: BigFloatRep.h
 * Synopsis:
 *                 Internal Representation BigFloat.
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

#ifndef _CORE_BIGFLOATREP_H_
#define _CORE_BIGFLOATREP_H_

#include <CGAL/CORE/BigRat.h>
#include <CGAL/CORE/CoreAux.h>
#include <CGAL/CORE/CoreDefs.h>
#include <CGAL/CORE/extLong.h>

namespace CORE {

//  forward reference
class BigFloat;

//  class BigFloatRep (internal representation for BigFloat)
class CGAL_CORE_EXPORT BigFloatRep : public RCRepImpl<BigFloatRep> {
public:
  static long chunkCeil(long bits);  //inline
  static long chunkFloor(long bits); //inline
  static long bits(long chunks);     //inline
  static BigInt chunkShift(const BigInt& x, long s); //inline
  static double lg10(BigInt x);      //inline
  static long floorlg10(BigInt x);   //inline

  /// exp2(e) returns 2^e : called by BigFloat::exp2(e)
  /** e can be negative */
  static BigFloatRep* exp2(int e);

  struct DecimalOutput;

  friend class BigFloat;

  BigInt        m;
  unsigned long err;
  long          exp;

public:
  //  constructors
  BigFloatRep(int=0);           //inline
  BigFloatRep(unsigned int);           //inline
  BigFloatRep(short);           //inline
  BigFloatRep(float);           //inline
  BigFloatRep(long);          //inline
  BigFloatRep(double);        //inline
  BigFloatRep(const BigInt& I, unsigned long u, long l); //inline
  BigFloatRep(const BigInt& I, long l); //inline
  BigFloatRep(const BigInt& I); //inline
  BigFloatRep(const char *);  //inline

  BigRat BigRatize() const;   //inline

  //  the destructor
  ~BigFloatRep(); //inline

  CORE_NEW(BigFloatRep)    // allocate the memory pool, unless
  CORE_DELETE(BigFloatRep) // memory pool feature is disabled.

  //  approximation
  void trunc(const BigInt&, const extLong&, const extLong&);
  void truncM(const BigFloatRep&, const extLong&, const extLong&);
  void approx(const BigFloatRep&, const extLong&, const extLong&);

  void div(const BigInt&, const BigInt&, const extLong&, const extLong&);
  void approx(const BigRat&, const extLong&, const extLong&); //inline

  //  error-normalization
  void eliminateTrailingZeroes(); //inline
  void normal();
  void bigNormal(BigInt&);

  //  arithmetics
public:
  void add(const BigFloatRep&, const BigFloatRep&);
  void sub(const BigFloatRep&, const BigFloatRep&);
  void mul(const BigFloatRep&, const BigFloatRep&);
  void div(const BigFloatRep&, const BigFloatRep&, const extLong&);
  void div2(const BigFloatRep&); ///< exact division by 2
  /// Converts a pair of BigFloatReps into one with error bounds
  void centerize(const BigFloatRep&, const BigFloatRep&);
private:
  //  squareroot
  //    arguments:      r = value whose square root we want
  //                    a = absolute precision of the desired result
  //                    init = initial approx. to the square root (for Newton)
  void sqrt(const BigInt& r, const extLong& a);
  /// sqrt(r,a,rr) -- compute sqrt(r) to absolute precision a,
  ///      starting from initial approximation of rr.
  void sqrt(const BigInt& r, const extLong& a, const BigFloat& init);
  void sqrt(const BigFloatRep& r, const extLong& a);
  /// sqrt(r,a,rr) -- compute sqrt(r) to absolute precision a,
  ///      starting from initial approximation of rr.
  void sqrt(const BigFloatRep& r, const extLong& a, const BigFloat& init);

  //  comparison
  int compareMExp(const BigFloatRep&) const;

  //  builtin functions
  extLong lMSB() const;      //inline
  extLong uMSB() const;      //inline
  extLong MSB() const;       //inline
  extLong flrLgErr() const;  //inline
  extLong clLgErr() const;   //inline

  bool    isZeroIn() const;  //inline
  int     signM() const;     //inline

  //  cast functions
  double toDouble() const;
  long toLong() const;
  BigInt toBigInt() const;

  //  conversion

  // toString() Joaquin Grech 31/5/2003
  std::string toString(long prec=get_static_defBigFloatOutputDigits(), bool sci=false) const;
  std::string round(std::string inRep, long& L10, unsigned int width) const;
  DecimalOutput toDecimal(unsigned int width=get_static_defBigFloatOutputDigits(),
                          bool Scientific=false) const;
  void fromString(const char *p, extLong  prec = getBigFloatInputDigits());

  void dump() const;  //inline
  long adjustE(long E, BigInt M, long e) const;

public:
  //  stream operators
  std::ostream& operator <<(std::ostream&) const; //inline
  std::istream& operator >>(std::istream&);
};//class BigFloatRep

////////////////////////////////////////////////////////////
//  Implementations
////////////////////////////////////////////////////////////

struct BigFloatRep::DecimalOutput {
  std::string rep;    // decimal output string
  int sign;           // 0, +1 or -1
  bool isScientific;  // false=positional notation
  int noSignificant;  // number of significant digits
  //   -1 means this information is not explicitly
  //   given, and must be determined from rep, etc.
  bool isExact;       //
  int errorCode;      // 0 = no error
                      // 1 = sign of number is unknown (e.g., mantissa
                      //  is smaller than error)

  DecimalOutput() : rep(""), sign(1), isScientific(false),
                    noSignificant(0), isExact(false), errorCode(0) {}
};//DecimalOutput

// constants used by BigFloatRep
//        NOTES:  CHUNK_BIT is the number of bits in each Chunk
//        Since LONG_BIT = 32 or 64, then CHUNK_BIT = 14 or 30.
//        We have:  0 <= err < 4 * 2^{CHUNK_BIT}

const long CHUNK_BIT = (long)(LONG_BIT / 2 - 2);         //  chunks
const long HALF_CHUNK_BIT = (CHUNK_BIT + 1) / 2;
const long DBL_MAX_CHUNK = (DBL_MAX_EXP - 1) / CHUNK_BIT + 1;
const double lgTenM = 3.321928094887362;

inline long BigFloatRep::chunkCeil(long bits) {
  if (bits > 0)
    return (bits - 1) / CHUNK_BIT + 1;
  else
    return - (- bits) / CHUNK_BIT;
}//chunkCeil

inline long BigFloatRep::chunkFloor(long bits) {
  if (bits >= 0)
    return bits / CHUNK_BIT;
  else
    return - (- bits - 1) / CHUNK_BIT - 1;
}//chunkFloor

// bits(c) returns the number of bits in c chunks:
inline long BigFloatRep::bits(long chunks) {
  return CHUNK_BIT * chunks;
}

inline BigInt BigFloatRep::chunkShift(const BigInt& x, long s) {
  if (!s || sign(x) == 0)
    return x;
  else if (s > 0)
    //  shift left
    if (sign(x) > 0)
      return x << static_cast<unsigned long>(bits(s));
    else //  x < 0
      return - ((-x) << static_cast<unsigned long>(bits(s)));
  else //  shift right
    if (sign(x) > 0)
      return x >> static_cast<unsigned long>(bits(-s));
    else //  x < 0
      return - ((-x) >> static_cast<unsigned long>(bits(-s)));
}//chunkShift

inline BigFloatRep* BigFloatRep::exp2(int e) {
  long ee;  // this is going to be the exponent
  if (e >= 0)
    ee = e/CHUNK_BIT;
  else
    ee = - ((-e + CHUNK_BIT -1)/CHUNK_BIT);

  int rem = e - (ee * CHUNK_BIT);     // Assert: 0 <= rem < CHUNK_BIT

  return new BigFloatRep((1<<rem), 0, ee);
  // Here, we assume CHUNK_BIT is less than int width
}

//  constructors
inline BigFloatRep::BigFloatRep(short n)
  : m(n), err(0), exp(0) {}

inline BigFloatRep::BigFloatRep(float n)
  : m(n), err(0), exp(0) {}

//  Chee (8/8/04) -- introduced constructor from int
inline BigFloatRep::BigFloatRep(int n)
  : m(n), err(0), exp(0) {}

inline BigFloatRep::BigFloatRep(unsigned int n)
  : m(n), err(0), exp(0) {}

//  Chee (8/8/04) -- introduced constructor from long
inline BigFloatRep::BigFloatRep(long n)
  : m(n), err(0), exp(0) {}

//  Chee (8/8/04) -- introduced constructor from double
/* This turns out to be an alternative implementation of the
 * original one in BigFloat.cpp!!
inline BigFloatRep::BigFloatRep(double d)
            : m(IntMantissa(d)), err(0), exp(0) {
     BigFloatRep * bfr = exp2(IntExponent(d));  // take care of the exponent
     m *= bfr->m;
     exp = bfr->exp;
}
*/

inline BigFloatRep::BigFloatRep(const BigInt& I, unsigned long er, long ex)
  : m(I), err(er), exp(ex) {}

inline BigFloatRep::BigFloatRep(const BigInt& I)
  : m(I), err(0), exp(0) {}

//Constructs the BigFloat representing I*2^{ex}.
//If ex >=0 then it is clear how to do it.
//Otherwise, let |ex| = CHUNK_BIT * q + r. Then
//I*2^{ex} = I*2^{CHUNK_BIT -r} 2^{-CHUNK_BIT * (q+1)}
inline BigFloatRep::BigFloatRep(const BigInt& I, long ex) {
  err=0;
  exp = chunkFloor(ex);
  if(ex >= 0)
    m = I<<(ex - bits(exp));
  else{//ex < 0
    exp = chunkFloor(abs(ex));
    m = I << (CHUNK_BIT - (-ex - bits(exp)));
    exp = -1*(1 + exp);
  }
}

inline BigFloatRep::BigFloatRep(const char *str) : m(0), err(0), exp(0) {
  fromString(str);
}

inline BigRat BigFloatRep::BigRatize() const {
  if (exp >= 0)
    return BigRat(chunkShift(m, exp), 1);
  else
    return BigRat(m, chunkShift(1, - exp));
}

//  the destructor
inline BigFloatRep::~BigFloatRep() {}

inline void BigFloatRep::approx(const BigRat& R, const extLong& r, const extLong& a) {
  div(numerator(R), denominator(R), r, a);
}

//  eliminate trailing zeroes
inline void BigFloatRep::eliminateTrailingZeroes() {
  // eliminate trailing 0's    -- IP 10/9/98
  /*if (err == 0 && m != 0) {
    while ((m & ((1 << CHUNK_BIT) - 1)) == 0) {
      m >>= CHUNK_BIT;
      exp++;
    }
  }*/
  // new code, much faster, Zilin Du (Nov, 2003)
  if (err == 0 && sign(m) != 0) {
    int r = getBinExpo(m) / CHUNK_BIT;
    m >>= (r * CHUNK_BIT);
    exp += r;
  }
}

//  bultin functions
inline extLong BigFloatRep::lMSB() const {
  if (!isZeroIn())
    return extLong(floorLg(abs(m) - err)) + bits(exp);
  else
    return extLong(CORE_negInfty);
}

/// uMSB() returns an upper bound on log_2(abs(*this)).
/** Returns -1 if (*this)=0.
 * Not well-defined if zero is in the interval.
 */
inline extLong BigFloatRep::uMSB() const {
  return extLong(floorLg(abs(m) + err)) + bits(exp);
}

inline extLong BigFloatRep::MSB() const {
  // Note : MSB is undefined if it's not exact.
  if (sign(m))          // sign(m) is non-zero
    return extLong(floorLg(m)) + bits(exp);
  else
    return extLong(CORE_negInfty);
}

inline extLong BigFloatRep::flrLgErr() const {
  if (err)
    return extLong(flrLg(err)) + bits(exp);
  else
    return extLong(CORE_negInfty);
}

inline extLong BigFloatRep::clLgErr() const {
  if (err)
    return extLong(clLg(err)) + bits(exp);
  else
    return extLong(CORE_negInfty);
}

// isZero() = true iff zero is inside the interval of BigFloat:
inline bool BigFloatRep::isZeroIn() const {
  if (err == 0){
    return (m == 0);        //Nov 6, 2002: bug fix!
  }
  long lm = bitLength(m);
  if (lm > CHUNK_BIT+2) {
    return false;   // since err < 4 * 2^{CHUNK_BIT}
  } else {
    return (abs(m) <= BigInt(err));
  }
}

inline int BigFloatRep::signM() const {
  return sign(m);
}

inline double BigFloatRep::lg10(BigInt x) {
  if (x == 0)
    return 0;

  BigInt t(abs(x));
  long l = -1;
  double d = 0;

  while (t > 0) {
    l++;
    d /= 10;
    d += ulongValue(t%10);
    t /= 10;
  }
  return std::log10(d) + l;
}

// this is a simpler form of lg10()
inline long BigFloatRep::floorlg10(BigInt x) {
  if (x == 0)
    return 0;
  BigInt t(abs(x));
  long l = -1;

  while (t > 0) {
    l++;
    t /= 10;
  }
  return l;
}

inline std::ostream& BigFloatRep::operator<<(std::ostream& o) const {
  bool sci = (o.flags() & std::ios::scientific) > 0;
  BigFloatRep::DecimalOutput r = toDecimal(o.precision(), sci);
  if (r.sign == -1)
    o << "-";
  o << r.rep;
  return o;
}

/* Returns a std::string with precision and format specified
   Works as cout << with the exception that if the output
   contains any error it returns a nullptr
   Joaquin Grech 31/5/03
   */
inline std::string BigFloatRep::toString(long prec, bool sci) const {
  BigFloatRep::DecimalOutput r = toDecimal(prec, sci);

  if (r.errorCode == 0) {
    if (r.sign < 0)
      return std::string("-")+r.rep;
    else
      return r.rep;
  }
  return nullptr;
}

inline void BigFloatRep::dump() const {
  std::cout << "---- BFRep: " << this << " ----" << std::endl;
  std::cout << "  BF value: ";
  this->operator<<(std::cout);
  std::cout <<  std::endl;
  std::cout << "  m = " << m << std::endl;
  std::cout << "  err = " << err << std::endl;
  std::cout << "  exp = " << exp << std::endl;
  std::cout << " -- End of BFRep " << this << " -- " << std::endl;
}



} //namespace CORE
#endif // _CORE_BIGFLOATREP_H_
