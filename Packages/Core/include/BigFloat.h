/******************************************************************
 * Core Library Version 1.5, August 2002
 * Copyright (c) 1995-2002 Exact Computation Project
 * 
 * File: BigFloat.h
 * Synopsis:
 *       An implementation of BigFloat numbers with error bounds.
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

#ifndef CORE_BIGFLOAT_H
#define CORE_BIGFLOAT_H

#ifdef DEBUG
#include <assert.h>
#endif

#include "CoreImpl.h"
#include "CoreAux.h"
#include "CoreDefs.h"
#include "extLong.h"
#include "BigInt.h"
#include "BigRat.h"
#include "MemoryPool.h"

CORE_BEGIN_NAMESPACE

//  forward reference

class BigFloat;

//  class BigFloatRep

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
  
  // exp2(e) returns 2^e : ?? Is it used?
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

//  class BigFloat
//  almost all functions are inline except the last four I/O functions.

class BigFloat
{
  BigFloatRep* rep;
  
  //  comparisons  
  int compare(const BigFloat&) const;
  BigFloat(BigFloatRep * r);
        
public:  
  //  constructors  
  BigFloat();
  BigFloat(int);
  BigFloat(long);
  BigFloat(double);
  BigFloat(const BigInt&, unsigned long = 0, long = 0);
  BigFloat(const char *);
  BigFloat(const BigRat & R,
        const extLong & r = defRelPrec,
        const extLong & a = defAbsPrec );
  BigFloat(const BigFloat&);
  
  BigInt m() const { return rep->m; }
  unsigned long err() const { return rep->err; }
  long exp() const { return rep->exp; }
  BigFloatRep * getRep() const { return rep; }

  static BigFloat exp2(int e) { // returns a BigFloat value 2^e
     return BigFloat(BigFloatRep::exp2(e));
  }


  // Chen: To resolve ambiguities between user-defined conversion operator
  // and user-defined constructors, we change the operator BigRat() const 
  // to the following function BigRatize() const:
  BigRat BigRatize() const;  // this function converts a BigFloat to a BigRat
 
  //  the destructor 
  ~BigFloat();
  
  //  assignment operator 
  BigFloat& operator= (const BigFloat&);
  
  //  approximation
  void approx(const BigInt& I, const extLong& r, const extLong& a);
  void approx(const BigFloat& B, const extLong& r, const extLong& a);
  void approx(const BigRat& R, const extLong& r, const extLong& a);
  
  //  unary minus
  BigFloat operator- () const;
  
  //  arithmetics  
  friend BigFloat operator+ (const BigFloat&, const BigFloat&);
  friend BigFloat operator- (const BigFloat&, const BigFloat&);
  friend BigFloat operator* (const BigFloat&, const BigFloat&);
  friend BigFloat operator/ (const BigFloat&, const BigFloat&);
  
  BigFloat div(const BigFloat& x, const extLong& r) const;
  
  //  squareroot
  BigFloat sqrt(const extLong&) const;
  BigFloat sqrt(const extLong&, const BigFloat& init) const;
                // init is initial approximation
  friend BigFloat sqrt(const BigFloat&);
  
  //  comparisons 
  friend bool operator== (const BigFloat&, const BigFloat&);
  friend bool operator!= (const BigFloat&, const BigFloat&);
  friend bool operator < (const BigFloat&, const BigFloat&);
  friend bool operator<= (const BigFloat&, const BigFloat&);
  friend bool operator > (const BigFloat&, const BigFloat&);
  friend bool operator>= (const BigFloat&, const BigFloat&);
  
  //  arithmetic and assignment opeartors
  BigFloat& operator +=(const BigFloat&);
  BigFloat& operator -=(const BigFloat&);
  BigFloat& operator *=(const BigFloat&);
  BigFloat& operator /=(const BigFloat&);

  //  cast operators [Sylvain, August 5, 2002 : remove the automatic conversion]
  double toDouble() const;
  float  toFloat() const;
  long   toLong() const;
  int    toInt() const;
  BigInt toBigInt() const;
  
  //  builtin function
  bool    isExact() const;	// true if err==0
  BigFloat& makeExact()      	// set err to 0 (return an exact BigFloat)
	{ rep->err =0; return *this;}   
  extLong lMSB() const;
  extLong uMSB() const;
  extLong MSB() const;
  extLong flrLgErr() const;
  extLong clLgErr() const;
  bool    isZeroIn() const;
  
  int 	  sign() const;		// This is just the sign of the mantissa!!
                                // This can be taken to be the sign of
                                // the BigFloat only if !(isZeroIn()).
  friend int sign(const BigFloat&); 	// again, this is only the sign
				// of the mantissa
  BigFloat abs() const;
  friend BigFloat abs(const BigFloat&);

  // helper function
  void dump() const;
  
  //  stream operators
  friend std::ostream& operator <<(std::ostream&, const BigFloat&);
  friend std::istream& operator >>(std::istream&, BigFloat&);
  
  void read_from_file(std::istream& in, long maxLength = 0);
  void read_from_file2(std::istream& in, long maxLength = 0);
  void write_to_file(std::ostream& out, int base = 10, int charPerLine = 80);
  void write_to_file2(std::ostream& out, int base = 10, int charPerLine = 80);
};

#ifdef CORE_ENABLE_INLINES
#include "BigFloat.inl"
#endif

CORE_END_NAMESPACE

#endif
