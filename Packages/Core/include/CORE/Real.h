/******************************************************************
 * Core Library Version 1.6, June 2003
 * Copyright (c) 1995-2002 Exact Computation Project
 * 
 * File: Real.h
 *
 * Synopsis: The Real class is a superclass for all the number 
 *           systems in the Core Library (int, long, float, double,
 *           BigInt, BigRat, BigFloat, etc)
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
 * $Id$
 *****************************************************************/

#ifndef CORE_REAL_H
#define CORE_REAL_H

#include "RealRep.h"

CORE_BEGIN_NAMESPACE

/// \class Real Real.h
/// \brief Real is a superclass for all the number systems 

class Real
{
public:
  RealRep* rep; ///< handle to the "real" representation

  /// \name Constructors and Destructor
  //@{
  /// default constructor
  Real() { rep = new RealLong(0); }
  /// copy constructor
  Real(const Real& x) { rep = x.rep; rep->refCount++; }
  
  /// constructor for <tt>int</tt>
  Real(int i) { rep = new RealLong(i); }
  /// constructor for <tt>unsigned int</tt>
  Real(unsigned int ui) {
    if (ui <= INT_MAX) // use RealLong, safe to convert ul to int.
      rep = new RealLong(static_cast<int>(ui));
    else // use BigInt since machine double has only 53-bit mantissa
      rep = new RealBigInt(ui);
  }
  
  /// constructor for <tt>long</tt>
  Real(long l) { rep = new RealLong(l); }
  /// constructor for <tt>unsigned long</tt>
  Real(unsigned long ul) {
    if (ul <= LONG_MAX) // use RealLong, safe to convert ul to long.
      rep = new RealLong(static_cast<long>(ul));
    else // use BigInt since machine double has only 53-bit mantissa
      rep = new RealBigInt(ul);
  }
  
  /// constructor for <tt>float</tt>
  Real(float f) { rep = new RealDouble(f); }
  /// constructor for <tt>double</tt>
  Real(double d) { rep = new RealDouble(d); }
  
  /// constructor for <tt>BigInt</tt>
  Real(const BigInt& I) { rep = new RealBigInt(I); }
  /// constructor for <tt>BigRat</tt>
  Real(const BigRat& R) { rep = new RealBigRat(R); }
  /// constructor for <tt>BigFloat</tt>
  Real(const BigFloat& B) { rep = new RealBigFloat(B); }

  /// constructor for <tt>const char *</tt>
  /** construct Real from a string representation \a s 
   * with precision \a prec */
  Real(const char *str, const extLong& prec = defInputDigits) : rep(NULL)
  { constructFromString(str, prec); }
  /// constructor for <tt>std::string</tt>
  Real(const std::string& s, const extLong& prec = defInputDigits) : rep(NULL)
  { constructFromString(s.c_str(), prec); }

  /// destructor
  ~Real() { if (--rep->refCount == 0) delete rep; }
  //@}

  /// \name Assignment Operators
  //@{ 
  /// assignment operator
  Real& operator=(const Real& x) {
  	if (this == &x)  return (*this); 
    if (--rep->refCount == 0) delete rep;
    rep = x.rep; rep->refCount++;
    return (*this);
  }
  //@}
  
  /// \name Compound Assignment Operators
  //@{ 
  /// operator+= 
  Real& operator+=(const Real& x) {
    Real t = *rep + x;
    if (--rep->refCount == 0) delete rep;
    rep = t.rep; rep->refCount++;
    return *this;
  }
  /// operator-= 
  Real& operator-=(const Real& x) {
    Real t = *rep - x;
    if (--rep->refCount == 0) delete rep;
    rep = t.rep; rep->refCount++;
    return *this;
  }
  /// operator*= 
  Real& operator*=(const Real& x) {
    Real t = *rep * x;
    if (--rep->refCount == 0) delete rep;
    rep = t.rep; rep->refCount++;
    return *this;
  }
  /// operator/= 
  Real& operator/=(const Real& x) {
  	Real t = rep->div(x, defRelPrec);
    if (--rep->refCount == 0) delete rep;
    rep = t.rep; rep->refCount++;
    return *this;
  }
  //@}

  /// \name Unary Minus, Increment and Decrement Operators
  //@{
  /// unary minus
  Real operator-() const { return -(*rep); }
  /// left increment operator (++i)
  Real& operator++() { *this += 1; return *this; }
  /// left decrement operator (--i)
  Real& operator--() { *this -= 1; return *this; }
  /// right increment operator (i++)
  Real operator++(int) { Real t = *this; *this += 1; return t; }
  /// right deccrement operator (i--)
  Real operator--(int) { Real t = *this; *this -= 1; return t; }
  //@}

  /// \name String Conversion Functions
  //@{
  /// set value from <tt>const char*</tt>
  void fromString(const char* s, const extLong& prec = defInputDigits)
  { *this = Real(s, prec); }
  /// convert to <tt>std::string</tt>
  /** give decimal string representation */
  std::string toString(long prec=defOutputDigits, bool sci=false) const
  { return rep->toString(prec, sci); }		  
  //@}

  /// \name Conversion Functions
  //@{
  /// convert to \c int
  int intValue() const { return (int) (*rep); }
  /// convert to \c long
  long longValue() const { return (long) (*rep); }
  /// convert to \c float
  float floatValue() const { return (float) (*rep); }
  /// convert to \c double
  double doubleValue() const { return (double) (*rep); }
  /// convert to \c BigInt
  BigInt BigIntValue() const { return rep->BigIntValue(); }
  /// convert to \c BigRat
  BigRat BigRatValue() const { return rep->BigRatValue(); }
  /// convert to \c BigFloat (approximate it first!)
  BigFloat BigFloatValue() const { return rep->BigFloatValue(); }    	
  //@}

  /// \name Aprroximation Function
  //@{
  /// approximation
  Real approx(const extLong& relPrec = defRelPrec,
        const extLong& absPrec = defAbsPrec) const
  { return rep->approx(relPrec, absPrec); }
  //@}

  /// \name Helper Functions
  //@{ 
  /// sign function
  int sign() const { return rep->sgn(); }
  /// isZero function
  bool isZero() { return sign() == 0; }
  /// return true if interval contains zero
  bool isZeroIn() const { return rep->isZeroIn(); }
  /// compare function
  int cmp(const Real& r) const { return (rep->operator-(r)).sign(); }
  /// absolute value function
  Real abs() const
  { Real x = (sign() >= 0) ? (*this) : -(*this); return x; }

  /// get mantissa of current approximate value
  BigInt getMantissa() const { return BigFloatValue().m(); }
  /// get exponent of current approximate value
  long getExponent() const { return BigFloatValue().exp(); }   
  /// return rep pointer:
  RealRep* getRep() const { return rep; }
  
  /// return true if error free otherwise return false; 
  bool  isExact() const { return rep->isExact(); }
  
  /// low bound of MSB
  extLong lMSB() const {
  	if (isExact()) { 
      return rep->mostSignificantBit;
    } else { // May 30, 2003:  Dirty Cast Bug Fix by Joaquin Grech
      //  -- this was causing crashes when called in Expr::withinKnownPrecision
      //  (necessitating -fno-strict-aliasing flag).
      BigFloat bf = rep->BigFloatValue();
      return bf.lMSB();
    }
  }
  /// upper bound of MSB
  extLong uMSB() const {
  	if (isExact()) {
      return rep->mostSignificantBit;
    } else {
      BigFloat bf = rep->BigFloatValue();
      return bf.uMSB();
    }
  }
  /// MSB - Most Significant Bit
  extLong MSB() const { return rep->mostSignificantBit; }
  
  /// floor of log_2 of Error
  extLong flrLgErr() const { return rep->flrLgErr(); }
  /// ceil of log_2 of Error
  extLong clLgErr() const { return rep->clLgErr(); }
  
  /// division with desired precision
  Real div(const Real& x, const extLong& r) const { return rep->div(x, r); }
  /// squareroot
  Real sqrt(const extLong& x) const { return rep->sqrt(x); }
  /// squareroot with initial approximation
  Real sqrt(const extLong& x, const BigFloat& A) const 
  { return rep->sqrt(x, A); }
  
  /// correspond to the variables "u25, l25, v2p, v2m, v5p, v5m" in Expr
  void ULV_E(extLong &up, extLong &lp, extLong &v2p, extLong &v2m,
	     extLong &v5p, extLong &v5m) const
  { rep->ULV_E(up, lp, v2p, v2m, v5p, v5m); }	     
  
  /// degree of polynomial P(x)
  unsigned long degree() const { return rep->degree(); }
  /// \f$ lg(|| P(X) ||_2) \f$
  unsigned long length() const { return rep->length(); }
  /// \f$ lg(|| P(X) ||_\infty) \f$
  unsigned long height() const { return rep->height(); }
  //@}

  /// return Real(0)
  static const Real& getZero();
  
private:
  void constructFromString(const char *str, const extLong& prec);
};//Class Real

#define CORE_REAL_ZERO Real::getZero()

/// I/O Stream operator<<
inline std::ostream& operator<< (std::ostream& o, const Real& r)
{ (r.rep)->operator <<(o); return o; }
/// I/O Stream operator>>
std::istream& operator >>(std::istream& i, Real& x);

/// addition
inline Real operator+ (const Real& x, const Real& y)
{ return x.rep->operator+ (y); }
/// substraction
inline Real operator- (const Real& x, const Real& y)
{ return x.rep->operator- (y); }
/// multiplication
inline Real operator* (const Real& x, const Real& y)
{ return x.rep->operator* (y); }
/// division
inline Real operator/ (const Real& x, const Real& y)
{ return x.rep->div(y, defRelPrec); }

/// operator ==
inline bool operator==(const Real& x, const Real& y)
{ return x.cmp(y) == 0; }
/// operator !=
inline bool operator!=(const Real& x, const Real& y)
{ return x.cmp(y) != 0; }
/// operator <
inline bool operator< (const Real& x, const Real& y)
{ return x.cmp(y) < 0; }
/// operator <=
inline bool operator<=(const Real& x, const Real& y)
{ return x.cmp(y) <= 0; }
/// operator >
inline bool operator> (const Real& x, const Real& y)
{ return x.cmp(y) > 0; }
/// operator >=
inline bool operator>=(const Real& x, const Real& y)
{ return x.cmp(y) >= 0; }

/// floor function
BigInt floor(const Real&, Real&);
/// power function
Real pow(const Real&, unsigned long);

/// return sign
inline int sign(const Real& r) { return r.sign(); }
/// is zero?
inline bool isZero(const Real& r) { return r.sign() == 0; }
/// compare
/** compare two Real \a r1 and \a r2, return
 * \retval -1 if r1 < r2,
 * \retval 0 if r1 = r2,
 * \retval 1 if r1 > r2. */
inline int cmp(const Real& r1, const Real& r2) { return r1.cmp(r2); }
/// absolute value
inline Real abs(const Real& x) { return x.abs(); }
/// absolute value (same as abs)
inline Real fabs(const Real& x) { return abs(x); }
/// floor
inline BigInt floor(const Real& r) { Real tmp; return floor(r, tmp); }
/// ceiling
inline BigInt ceil(const Real& r) { return -floor(-r); }
/// power
inline Real power(const Real& r, unsigned long p) { return pow(r, p); }
/// square root
inline Real sqrt(const Real& x) { return x.sqrt(defAbsPrec); }

#ifdef CORE_ENABLE_INLINES
  #include "Real.inl"
#endif

CORE_END_NAMESPACE
#endif // CORE_REAL_H
