/******************************************************************
 * Core Library Version 1.6, June 2003
 * Copyright (c) 1995-2003 Exact Computation Project
 * 
 * File: BigFloat.h
 *
 * Synopsis: An implementation of BigFloat numbers with error bounds.
 *
 * Written by 
 *       Chee Yap <yap@cs.nyu.edu>
 *       Chen Li  <chenli@cs.nyu.edu>
 *       Zilin Du <zilin@cs.nyu.edu>
 *
 * WWW URL: http://cs.nyu.edu/exact/
 * Email: exact@cs.nyu.edu
 *
 * $Id$
 *****************************************************************/

#ifndef CORE_BIGFLOAT_H
#define CORE_BIGFLOAT_H

#include "BigFloatRep.h"

CORE_BEGIN_NAMESPACE

/// \class BigFloat BigFloat.h
/// \brief BigFloat is a class of Float-Point number with error bounds.
class BigFloat
{
public:  
  BigFloatRep* rep; ///< handle to the "real" representation

  /// \name Constructors and Destructor
  //@{
  /// default constructor
  BigFloat() { rep = new BigFloatRep(); rep->refCount++; }
  /// copy constructor
  BigFloat(const BigFloat& x) { rep = x.rep; rep->refCount++; }
  
  /// constructor for <tt>int</tt>
  BigFloat(int i) { rep = new BigFloatRep(i); rep->refCount++; }
  /// constructor for <tt>long</tt>
  BigFloat(long l) { rep = new BigFloatRep(l); rep->refCount++; }

  /// constructor for <tt>double</tt>
  BigFloat(double d) { rep = new BigFloatRep(d); rep->refCount++; }

  /// constructor for <tt>const char* </tt>(default base = 10)
  BigFloat(const char *s) { rep = new BigFloatRep(s); rep->refCount++; }  
  /// constructor for <tt>std::string</tt>(default base = 10)
  BigFloat(const std::string& s) 
  { rep = new BigFloatRep(s.c_str()); rep->refCount++; }  

  /// constructor for <tt>BigInt</tt>
  BigFloat(const BigInt& I, unsigned long u = 0, long l = 0)
  { rep = new BigFloatRep(I, u, l); rep->refCount++; }
  /// constructor for <tt>BigRat</tt>
  BigFloat(const BigRat& R, const extLong& r = defRelPrec,
           const extLong& a = defAbsPrec)
  { rep = new BigFloatRep(); rep->refCount++; rep->approx(R, r, a); }
  
  /// constructor for <tt>BigFloatRep</tt>
  BigFloat(BigFloatRep * r) { this->rep = r; rep->refCount++; }
  
  /// destructor 
  ~BigFloat() { if (--rep->refCount == 0) delete rep; }
  //@}

  /// \name Assignment Operator
  //@{
  /// assignment operator 
  BigFloat& operator= (const BigFloat& x) {
    if (this == &x) return *this; 
    if (--rep->refCount == 0) delete rep;
    rep = x.rep; rep->refCount++;
    return *this;
  }
  //@}
  
  /// \name Compound Assignment Operators
  //@{
  /// operator+=
  BigFloat& operator+= (const BigFloat& x) {
    BigFloat z; z.rep->add(*rep, *x.rep);
    if (--rep->refCount == 0) delete rep;
    rep = z.rep; rep->refCount++;
    return *this;
  }
  /// operator-=  
  BigFloat& operator-= (const BigFloat& x) {
    BigFloat z; z.rep->sub(*rep, *x.rep);
    if (--rep->refCount == 0) delete rep;
    rep = z.rep; rep->refCount++;
    return *this;
  }
  /// operator*=  
  BigFloat& operator*= (const BigFloat& x) {
    BigFloat z; z.rep->mul(*rep, *x.rep);
    if (--rep->refCount == 0) delete rep;
    rep = z.rep; rep->refCount++;
    return *this;
  }
  /// operator/=  
  BigFloat& operator/= (const BigFloat& x) {
    BigFloat z; z.rep->div(*rep, *x.rep, defBFdivRelPrec);
    if (--rep->refCount == 0) delete rep;
    rep = z.rep; rep->refCount++;
    return *this;
  }
  //@}
  
  /// \name Unary Minus Operator
  //@{
  /// unary minus
  BigFloat operator- () const { return BigFloat(-rep->m, rep->err, rep->exp); }
  //@}
  
  /// \name String Conversion Functions
  //@{
  /// set value from <tt>const char*</tt> (base = 10)
  void fromString(const char* s, const extLong& p=defBigFloatInputDigits)
  { rep->fromString(s, p); }
  /// convert to <tt>std::string</tt> (base = 10)
  std::string toString(long prec=defBigFloatOutputDigits, bool sci=false) const 
  { return rep->toString(prec, sci); }
  //@}

  /// \name Conversion Functions
  //@{
  /// return int value
  int intValue() const { return (int)rep->toLong(); }
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
  float floatValue() const { return (float)rep->toDouble(); } 
  /// return double value
  double doubleValue() const { return rep->toDouble(); } 
  /// return BigInt value
  BigInt BigIntValue() const { return rep->toBigInt(); }
  /// return BigRat value
  BigRat BigRatValue() const { return rep->BigRatize(); }
  //@}
  
  /// \name File I/O Functions
  //@{
  /// read from file
  void readFromFile(std::istream& in, long maxLength = 0);
  /// write to file
  void writeToFile(std::ostream& in, int base=10, int charsPerLine=80) const;
  //@}
   
  /// \name Helper Functions
  //@{
  /// sign function
  /** \note This is only the sign of the mantissa, it can be taken to be 
      the sign of the BigFloat only if !(isZeroIn()). */
  int sign() const { return rep->signM(); }
  /// check whether contains zero
  /** \return true if contains zero, otherwise false */
  bool isZeroIn() const { return rep->isZeroIn(); }  
  /// absolute value function
  BigFloat abs() const 
  { return (sign()>=0) ? BigFloat(*this):BigFloat(-rep->m,rep->err,rep->exp);}
  ///  comparison function 
  int cmp(const BigFloat& x) const { return rep->compareMExp(*x.rep); }

  /// get mantissa    
  const BigInt& m() const { return rep->m; }
  /// get error bits
  unsigned long err() const { return rep->err; }
  /// get exponent
  long exp() const { return rep->exp; }
  /// return rep
  BigFloatRep* getRep() const { return rep; }
  
  /// check whether err == 0
  /** \return true if err == 0, otherwise false */
  bool isExact() const { return rep->err == 0; }
  /// set err to 0 
  /** \return an exact BigFloat, see Tutorial for why this is useful! */
  BigFloat& makeExact() { rep->err =0; return *this;}
  /// set err to 1
  /** \return an inexact BigFloat, see Tutorial for why this is useful! */
  BigFloat& makeInexact() { rep->err =1; return *this;}
  
  /// return lower bound of Most Significant Bit
  extLong lMSB() const { return rep->lMSB(); }
  /// return upper bound of Most Significant Bit
  extLong uMSB() const { return rep->uMSB(); }
  /// return Most Significant Bit
  extLong MSB() const { return rep->MSB(); }
  
  /// floor of Lg(err)
  extLong flrLgErr() const { return rep->flrLgErr(); }
  /// ceil of Lg(err)
  extLong clLgErr() const { return rep->clLgErr(); }  
  
  /// division with relative precsion <tt>r</tt>
  BigFloat div(const BigFloat& x, const extLong& r) const
  { BigFloat y; y.rep->div(*rep, *x.rep, r); return y; }
  /// exact division by 2
  BigFloat div2() const
  { BigFloat y; y.rep->div2(*rep); return y; }
  
  /// squareroot
  BigFloat sqrt(const extLong& a) const
  { BigFloat x;  x.rep->sqrt(*rep, a); return x; }
  /// squareroot with initial approximation <tt>init</tt>
  BigFloat sqrt(const extLong& a, const BigFloat& init) const
  { BigFloat x;  x.rep->sqrt(*rep, a, init); return x; }
  //@}

  /// \name Utility Functions
  //@{
  /// approximate BigInt number
  void approx(const BigInt& I, const extLong& r, const extLong& a) {
    if ((rep->refCount) > 1) {//  *rep is shared
      --rep->refCount;
      rep = new BigFloatRep();
      rep->refCount++;
    }
    rep->trunc(I, r, a);
  }
  /// approximate BigFloat number
  void approx(const BigFloat& B, const extLong& r, const extLong& a) {
    if ((rep->refCount) > 1) {//  *rep is shared
      --rep->refCount;
      rep = new BigFloatRep();
      rep->refCount++;
    }
    rep->approx(*B.rep, r, a);
  }
  /// approximate BigRat number
  void approx(const BigRat& R, const extLong& r, const extLong& a) {
    if ((rep->refCount) > 1) {//  *rep is shared
      --rep->refCount;
      rep = new BigFloatRep();
      rep->refCount++;
    }
    rep->approx(R, r, a);
  }
  /// dump internal data
  void dump() const { rep->dump(); }
  //@}

  /// returns a BigFloat of value \f$ 2^e \f$
  static BigFloat exp2(int e) { return BigFloat(BigFloatRep::exp2(e)); }
};

/// IO stream operator<<
inline std::ostream& operator<< (std::ostream& o, const BigFloat& x)
{ x.rep->operator<<(o); return o; }
/// IO stream operator>>
inline std::istream& operator>> (std::istream& i, BigFloat& x)
{ x.rep->operator>>(i); return i; }

/// operator+
inline BigFloat operator+ (const BigFloat& x, const BigFloat& y)
{ BigFloat z; z.rep->add(*x.rep, *y.rep); return z; }
/// operator-
inline BigFloat operator- (const BigFloat& x, const BigFloat& y)
{ BigFloat z; z.rep->sub(*x.rep, *y.rep); return z; }
/// operator*
inline BigFloat operator* (const BigFloat& x, const BigFloat& y)
{ BigFloat z; z.rep->mul(*x.rep, *y.rep); return z; }
/// operator/
inline BigFloat operator/ (const BigFloat& x, const BigFloat& y)
{ BigFloat z; z.rep->div(*x.rep, *y.rep, defBFdivRelPrec); return z; }

/// operator==
inline bool operator== (const BigFloat& x, const BigFloat& y)
{ return x.cmp(y) == 0; }
/// operator!=
inline bool operator!= (const BigFloat& x, const BigFloat& y)
{ return x.cmp(y) != 0; }
/// operator>=
inline bool operator>= (const BigFloat& x, const BigFloat& y)
{ return x.cmp(y) >= 0; }
/// operator>
inline bool operator> (const BigFloat& x, const BigFloat& y)
{ return x.cmp(y) > 0; }
/// operator<=
inline bool operator<= (const BigFloat& x, const BigFloat& y)
{ return x.cmp(y) <= 0; }
/// operator<
inline bool operator< (const BigFloat& x, const BigFloat& y)
{ return x.cmp(y) < 0; }

/// sign
inline int sign(const BigFloat& x) { return x.sign(); }
/// abs
inline BigFloat abs(const BigFloat& x) { return x.abs(); }
/// cmp
inline int cmp(const BigFloat& x, const BigFloat& y) { return x.cmp(y); }
/// pow
BigFloat pow(const BigFloat&, unsigned long);
/// power
inline BigFloat power(const BigFloat& x, unsigned long p) { return pow(x, p); }
/// sqrt to defAbsPrec:
inline BigFloat sqrt(const BigFloat& x) { return x.sqrt(defBFsqrtAbsPrec); }

/// convert an BigFloat Interval to a BigFloat with error bits
inline BigFloat centerize(const BigFloat& a, const BigFloat& b)
{ BigFloat z; z.rep->centerize(*a.rep, *b.rep); return z;}

CORE_END_NAMESPACE
#endif // CORE_BIGFLOATT_H
