/******************************************************************
 * Core Library Version 1.6, June 2003
 * Copyright (c) 1995-2003 Exact Computation Project
 * 
 * File: BigFloatRep.h
 *
 * Synopsis: Internal Representation BigFloat.
 *
 * Written by 
 *       Chee Yap <yap@cs.nyu.edu>
 *       Zilin Du <zilin@cs.nyu.edu>
 *
 * WWW URL: http://cs.nyu.edu/exact/
 * Email: exact@cs.nyu.edu
 *
 * $Id$
 *****************************************************************/

#ifndef CORE_BIGFLOATREP_H
#define CORE_BIGFLOATREP_H

#ifdef DEBUG
#include <assert.h>
#endif

#include <CORE/CoreImpl.h>
#include <CORE/CoreAux.h>
#include <CORE/CoreDefs.h>
#include <CORE/extLong.h>
#include <CORE/BigInt.h>
#include <CORE/BigRat.h>
#include <CORE/MemoryPool.h>

CORE_BEGIN_NAMESPACE

//  forward reference
class BigFloat;

//  class BigFloatRep (internal representation for BigFloat)
class BigFloatRep
{
private:
  static long chunkCeil(long bits);  //inline
  static long chunkFloor(long bits); //inline
  static long bits(long chunks); //inline
  static BigInt chunkShift(const BigInt& x, long s); //inline
  static double lg10(BigInt x); //inline
  static long floorlg10(BigInt x); //inline
  static void error(const char*); //inline
  
  /// exp2(e) returns 2^e : called by BigFloat::exp2(e)
  /** e can be negative */
  static BigFloatRep* exp2(int e);
  
  struct DecimalOutput;

  friend class BigFloat;
  
  BigInt        m;
  unsigned long err;
  long          exp;
  
  unsigned      refCount;  

public:  
  //  constructors
  BigFloatRep(double);
  BigFloatRep(const BigInt& I = 0, unsigned long u = 0, long l = 0); //inline
  BigFloatRep(const char *);  //inline
  
  BigRat BigRatize() const;   //inline
  
  //  the destructor
  ~BigFloatRep(); //inline

  CORE_MEMORY(BigFloatRep)
  
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
  void div2(const BigFloatRep&);	// exact division by 2
  void centerize(const BigFloatRep&, const BigFloatRep&);
private:  
  //  squareroot
  //    arguments:      r = value whose square root we want
  //                    a = absolute precision of the desired result
  //                    init = initial approx. to the square root (for Newton)
  void sqrt(const BigInt& r, const extLong& a);
  void sqrt(const BigInt& r, const extLong& a, const BigFloat& init);
  void sqrt(const BigFloatRep& r, const extLong& a);
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
  std::string toString(long prec=defBigFloatOutputDigits, bool sci=false) const;
  std::string round(std::string inRep, unsigned int width) const;
  DecimalOutput toDecimal(unsigned int width=defBigFloatOutputDigits, 
                          bool Scientific=false) const;
  void fromString(const char *p, const extLong & prec = defBigFloatInputDigits);
  
  void dump() const;  //inline
  long adjustE(long E, BigInt M, long e) const;

public:
  //  stream operators
  std::ostream& operator <<(std::ostream&) const; //inline
  std::istream& operator >>(std::istream&);
};

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
                    noSignificant(0), isExact(false), errorCode(0) 
  {}
};

// constants used by BigFloatRep
//	NOTES:  CHUNK_BIT is the number of bits in each Chunk
//	Since LONG_BIT = 32 or 64, then CHUNK_BIT = 14 or 30.  
//	We have:  0 <= err < 4 * 2^{CHUNK_BIT} 

const long CHUNK_BIT = (long)(LONG_BIT / 2 - 2); 	//  chunks
const long HALF_CHUNK_BIT = (CHUNK_BIT + 1) / 2;
const long DBL_MAX_CHUNK = (DBL_MAX_EXP - 1) / CHUNK_BIT + 1;
const double lgTenM = 3.321928094887362;

#ifdef CORE_ENABLE_INLINES
#include <CORE/BigFloat.inl>
#endif

CORE_END_NAMESPACE

#endif
