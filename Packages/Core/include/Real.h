/******************************************************************
 * Core Library Version 1.5, August 2002
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

#include "CoreImpl.h"
#include "CoreAux.h"
#include "BigFloat.h"
#include "MemoryPool.h"

CORE_BEGIN_NAMESPACE

//  forward reference
class RealRep;

/**  
 * \class Real Real.h
 * \brief Real is a superclass for all the number systems 
*/

class Real
{
public:
  /// \name Constructors and Destructor
  //@{
  /// default constructor
  Real();
  /// constructor for \c int
  Real(int);
  /// constructor for <tt>unsigned int</tt>
  Real(unsigned int);
  /// constructor for \c long
  Real(long);
  /// constructor for <tt>unsigned long</tt>
  Real(unsigned long);
  /// constructor for \c float
  Real(float);
  /// constructor for \c double
  Real(double);
  /// constructor for \c BigInt
  Real(BigInt);
  /// constructor for \c BigRat
  Real(BigRat);
  /// constructor for \c BigFloat
  Real(BigFloat);
  /// constructor for \c string
  /** construct Real from a string representation \a s 
   * with precision \a prec */
  Real(const char *str, const extLong& prec = defInputDigits);
  /// copy constructor
  Real(const Real&);
  /// destructor
  ~Real();
  //@}

  /// \name Aprroximation Function
  //@{
  // approximation
  Real approx(const extLong& relPrec = defRelPrec,
        const extLong& absPrec = defAbsPrec) const;
  //@}

  /// \name Help Functions
  //@{ 
  /// get current approximate value
  BigFloat getBigFloat() const;  	
  /// get exponent of current approximate value
  long getExponent() const; 
  /// get mantissa of current approximate value
  BigInt getMantissa() const;
  //@}

  /// \name Conversion Functions
  //@{
  /// convert to \c int
  int toInt() const ;
  /// convert to \c long
  long toLong() const ;
  /// convert to \c float
  float toFloat() const;
  /// convert to \c double
  double toDouble() const;
  /// convert to \c string
  /** give decimal string representation */
  const char *toString() const;
  /// convert to \c BigRat
  BigRat toBigRat() const;
  //@}

  /// \name Assignment Operators
  //@{ 
  /// assignment operator
  Real& operator=(const Real&);
  /// += operator
  Real& operator+=(const Real&);
  /// -= operator
  Real& operator-=(const Real&);
  /// *= operator
  Real& operator*=(const Real&);
  /// /= operator
  Real& operator/=(const Real&);
  //@}

  /// \name Increment, Decrement, and Unary Minus Operators
  //@{
  /// left increment operator (++i)
  Real& operator++();
  /// right increment operator (i++)
  Real operator++(int);
  /// left decrement operator (--i)
  Real& operator--();
  /// right deccrement operator (i--)
  Real operator--(int);
  /// unary minus
  Real operator-() const;
  //@}
  
  /// \name Arithematic Operators
  //@{
  /// addition
  friend Real operator+(const Real&, const Real&);
  /// substraction
  friend Real operator-(const Real&, const Real&);
  /// multiplication
  friend Real operator*(const Real&, const Real&);
  /// division
  friend Real operator/(const Real&, const Real&);
  /// square root
  friend Real sqrt(const Real&);
  //@}

  /// \name Comparison Operators
  //@{
  /// operator ==
  friend bool operator==(const Real&, const Real&);
  /// operator !=
  friend bool operator!=(const Real&, const Real&);
  /// operator <
  friend bool operator< (const Real&, const Real&);
  /// operator <=
  friend bool operator<=(const Real&, const Real&);
  /// operator <
  friend bool operator> (const Real&, const Real&);
  /// operator >=
  friend bool operator>=(const Real&, const Real&);
  //@}
  
  /// \name Builtin Functions
  //@{
  /// sign function
  friend int sign(const Real&);
  /// isZero function
  friend bool isZero(const Real&);
  /// compare function
  /** compare two Real \a e1 and \a e2, return
   * \retval -1 if e1 < e2,
   * \retval 0 if e1 = e2,
   * \retval 1 if e1 > e2. */
  friend int compare(const Real& e1, const Real& e2);
  /// floor function
  friend BigInt floor(const Real&);
  /// ceil function
  friend BigInt ceil(const Real&);
  /// power function
  friend Real pow(const Real&, unsigned long);
  /// power function (same as pow())
  friend Real power(const Real&, unsigned long n);
  /// absolute value function
  friend Real abs(const Real&);
  /// absolute value function (same as abs())
  friend Real fabs(const Real&);
  //@}

  /// \name I/O Stream
  //@{
  /// write to ostream 
  friend std::ostream& operator<<(std::ostream&, const Real&);
  /// read from istream
  friend std::istream& operator>>(std::istream&, Real&);
  //@}

public:
  /// \name Deprecated Functions
  //@{
  /// sign function
  int sign() const;
  /// return the pointer:
  RealRep* get_rep() const { return rep; }
  /// division with desired precision
  Real div(const Real&, const extLong&) const;
  /// squareroot
  Real sqrt(const extLong&) const; 
  /// squareroot with initial approximation
  Real sqrt(const extLong&, const BigFloat &) const; 

  /// return true if error free otherwise return false; 
  bool    isExact() const;
  /// low bound of MSB
  extLong lMSB() const;
  /// upper bound of MSB
  extLong uMSB() const;
  /// MSB - Most Significant Bit
  extLong MSB() const;
  
  /// correspond to the variables "u25, l25, v2p, v2m, v5p, v5m" in Expr
  void ULV_E(extLong &up, extLong &lp, extLong &v2p, extLong &v2m,
	     extLong &v5p, extLong &v5m) const;
  /// floor of log_2 of Error
  extLong flrLgErr() const;
  /// ceil of log_2 of Error
  extLong clLgErr() const;
  /// return true if interval contains zero
  bool    isZeroIn() const;
  /// degree of polynomial P(x)
  unsigned long degree() const;
  /// || P(X) ||_2
  unsigned long length() const;
  /// || P(X) ||_\infty
  unsigned long height() const;
  //@}
  /// return Real(0)
  static const Real& getZero();

protected:
  RealRep* rep; ///< handle to the "real" representation
};//Class Real


template <class T>
class Realbase_for;

typedef Realbase_for<long>     RealLong;
typedef Realbase_for<double>   RealDouble;
typedef Realbase_for<BigInt>   RealBigInt;
typedef Realbase_for<BigRat>   RealBigRat;
typedef Realbase_for<BigFloat> RealBigFloat;



/// \class RealRep 
/// \brief The internal representation of Real
class RealRep
{
public:
  /// reference counter
  unsigned refCount;
  
  /// most Significant Bit
  /** DEFINITION: MSB(0)=-\infty.  When E is not 0, and real, we define MSB(E)
      to be the natural number m such that 2^{m} <= |E| < 2^{m+1}. Hence, 
      MSB(E) is equal to floor(log_2(|E|)).  Intuitively, MSB is the position
      of the most significant bit in a binary rational representation of |E|.
      Thus, the bit before the binary point is considered to be position 0, 
      the the bit after the binary point is position -1.  Thus,
               ... b_2 b_1 b_0.b_{-1} b_{-2} ...
      E.g., MSB(1) = 0,
            MSB(1/2) = MSB((0.1)_2) = -1
	    MSB(1/4) = MSB((0.01)_2) = -2.
	    MSB(2) = MSB(3) = 1, MSB(4) = 2.
      Hence, if E is a non-zero integer, MSB(E) is equal to bitlength(|E|)-1.
      We also need upper and lower bounds on MSB(E).  This is defined to be
      any numbers lMSB(E) and uMSB(E) such that  lMSB(E) <= MSB(E) <= uMSB(E).
      THIS implies the following inequality:
            2^{lMSB(E)} <= |E| < 2^{1+uMSB(E)}.
      When E is an interval (e.g., BigFloat with non-zero error), then MSB(E)
      is not defined, but uMSB(E) and lMSB(E) is defined as follows:
      Assume E = [a, b].
         If 0 <  a <= b, then    lMSB(E) <= MSB(a) <= MSB(b) <= uMSB(E)
	 If a <= b <  0, then    lMSB(E) <= MSB(b) <= MSB(a) <= uMSB(E)
         If a <= 0 <= b, then    lMSB(E) = MSB(0) = -\infty
                                 uMSB(E) >=  max( MSB(a), MSB(b) ) */
  extLong mostSignificantBit;
  
  /// \name Constructor and Destructor
  //@{
  /// Constructor
  RealRep() : refCount(1) {}
  /// Destructor
  virtual ~RealRep() {}
  //@}
 
  /// \name Cast Operators
  //@{
  virtual operator double() const = 0;
  virtual operator float() const = 0;
  virtual operator long() const = 0;
  virtual operator int() const = 0;
  virtual BigFloat getBigFloat() const = 0;
  virtual BigRat toBigRat() const = 0;
  //@}
  
  /// \name Approximation
  //@{
  virtual Real approx(const extLong&, const extLong&) const = 0;
  //@}
 
  /// \name unary minus
  //@{
  virtual Real operator -() const = 0;
  //@}
  
  /// \name addition
  //@{
  virtual Real operator +(const Real&) const = 0;
  virtual Real addLong(const RealLong&) const = 0;
  virtual Real addDouble(const RealDouble&) const = 0;
  virtual Real addBigInt(const RealBigInt&) const = 0;
  virtual Real addBigFloat(const RealBigFloat&) const = 0;
  virtual Real addBigRat(const RealBigRat&) const = 0;
  //@}

  /// \name subtraction
  //@{
  virtual Real operator -(const Real&) const = 0;
  virtual Real subLong(const RealLong&) const = 0;
  virtual Real subDouble(const RealDouble&) const = 0;
  virtual Real subBigInt(const RealBigInt&) const = 0;
  virtual Real subBigFloat(const RealBigFloat&) const = 0;
  virtual Real subBigRat(const RealBigRat&) const = 0;
  //@}

  /// \name multiplication
  //@{
  virtual Real operator *(const Real&) const = 0;
  virtual Real mulLong(const RealLong&) const = 0;
  virtual Real mulDouble(const RealDouble&) const = 0;
  virtual Real mulBigInt(const RealBigInt&) const = 0;
  virtual Real mulBigFloat(const RealBigFloat&) const = 0;
  virtual Real mulBigRat(const RealBigRat&) const = 0;
  //@}  
  
  /// \name division
  //@{
  // virtual Real operator /(const Real&) const;
  virtual Real div(const Real&, const extLong&) const = 0;
  virtual Real divLong(const RealLong&, const extLong&) const = 0;
  virtual Real divDouble(const RealDouble&, const extLong&) const = 0;
  virtual Real divBigInt(const RealBigInt&, const extLong&) const = 0;
  virtual Real divBigFloat(const RealBigFloat&, const extLong&) const = 0;
  virtual Real divBigRat(const RealBigRat&, const extLong&) const = 0;
  //@}

  /// \name squareroot
  //@{
  virtual Real sqrt(const extLong&) const = 0;
  // sqrt with initial approximation
  virtual Real sqrt(const extLong&, const BigFloat &) const = 0; 
  //@}

  /// \name equality
  //@{
  virtual bool operator ==(const Real&) const = 0;
  virtual bool eqlLong(const RealLong&) const = 0;
  virtual bool eqlDouble(const RealDouble&) const = 0;
  virtual bool eqlBigInt(const RealBigInt&) const = 0;
  virtual bool eqlBigFloat(const RealBigFloat&) const = 0;
  virtual bool eqlBigRat(const RealBigRat&) const = 0;
  //@}
  
  /// \name smaller than
  //@{
  virtual bool operator <(const Real&) const = 0;
  virtual bool grtLong(const RealLong&) const = 0;
  virtual bool grtDouble(const RealDouble&) const = 0;
  virtual bool grtBigInt(const RealBigInt&) const = 0;
  virtual bool grtBigFloat(const RealBigFloat&) const = 0;
  virtual bool grtBigRat(const RealBigRat&) const = 0;
  //@}
  
  //  builtin functions
  //@{
  virtual bool    isExact() const = 0;
  virtual void ULV_E(extLong &up, extLong &lp, extLong &v2p, extLong &v2m,
		     extLong &v5p, extLong &v5m) const = 0;
  virtual extLong flrLgErr() const = 0;
  virtual extLong clLgErr() const = 0;
  virtual bool    isZeroIn() const = 0;
  virtual unsigned long degree() const = 0;
  virtual unsigned long length() const = 0;
  virtual unsigned long height() const = 0;
  virtual int sgn() const = 0;
  //@}
  //  I/O Stream
  //@{
  virtual std::ostream& operator <<(std::ostream&) const = 0;
  //@}
}; //class RealRep


template <class T>
class Realbase_for : public RealRep
{
public:
  CORE_MEMORY(Realbase_for)

  /// Kernel
  T ker;
  
  /// \name Constructor and Destructor
  //@{
  /// Constructor
  Realbase_for(const T&);
  /// Destructor
  ~Realbase_for() {}
  //@}
  
  /// Access to ker
  const T& get_ker() const { return ker; }

  /// \name cast operators
  //@{
  operator double() const;
  operator float() const;
  operator long() const;
  operator int() const;
  BigFloat getBigFloat() const;  
  BigRat toBigRat() const;  
  //@}

  /// \name approximation
  //@{
  Real approx(const extLong&, const extLong&) const;
  //@}
  
  /// \name unary minus
  //@{
  Real operator -() const;
  //@}
  
  /// \name addition
  //@{
  Real operator +(const Real&) const;
  Real addLong(const RealLong&) const;
  Real addDouble(const RealDouble&) const;
  Real addBigInt(const RealBigInt&) const;
  Real addBigFloat(const RealBigFloat&) const;
  Real addBigRat(const RealBigRat&) const;
  //@}
  
  /// \name subtraction
  //@{
  Real operator -(const Real&) const;
  Real subLong(const RealLong&) const;
  Real subDouble(const RealDouble&) const;
  Real subBigInt(const RealBigInt&) const;
  Real subBigFloat(const RealBigFloat&) const;
  Real subBigRat(const RealBigRat&) const;
  //@}
  
  /// \name multiplication
  //@{
  Real operator *(const Real&) const;
  Real mulLong(const RealLong&) const;
  Real mulDouble(const RealDouble&) const;
  Real mulBigInt(const RealBigInt&) const;
  Real mulBigFloat(const RealBigFloat&) const;
  Real mulBigRat(const RealBigRat&) const;
  //@}

  /// \name division
  //@{
  Real div(const Real&, const extLong&) const;
  Real divLong(const RealLong&, const extLong&) const;
  Real divDouble(const RealDouble&, const extLong&) const;
  Real divBigInt(const RealBigInt&, const extLong&) const;
  Real divBigFloat(const RealBigFloat&, const extLong&) const;
  Real divBigRat(const RealBigRat&, const extLong&) const;
  //@}
  
  /// \name squareroot
  //@{
  Real sqrt(const extLong&) const;
  Real sqrt(const extLong&,  const BigFloat&) const;
  //@}
  
  /// \name equality
  //@{
  bool operator ==(const Real&) const;
  bool eqlLong(const RealLong&) const;
  bool eqlDouble(const RealDouble&) const;
  bool eqlBigInt(const RealBigInt&) const;
  bool eqlBigFloat(const RealBigFloat&) const;
  bool eqlBigRat(const RealBigRat&) const;
  //@}
  
  /// \name smaller-than
  //@{
  bool operator <(const Real&) const;
  bool grtLong(const RealLong&) const;
  bool grtDouble(const RealDouble&) const;
  bool grtBigInt(const RealBigInt&) const;
  bool grtBigFloat(const RealBigFloat&) const;
  bool grtBigRat(const RealBigRat&) const;
  //@}
  
  /// \name builtin functions
  //@{
  bool    isExact() const;
  void ULV_E(extLong &up, extLong &lp, extLong &v2p, extLong &v2m,
	     extLong &v5p, extLong &v5m) const;
  extLong flrLgErr() const;
  extLong clLgErr() const;
  bool    isZeroIn() const; 
  unsigned long degree() const;
  unsigned long length() const;
  unsigned long height() const;
  int sgn() const;
  //@}

  /// \name I/O stream
  //@{
  std::ostream& operator <<(std::ostream&) const;
  //@}
};

BigInt floor(const Real&, Real&);
Real pow(const Real&, unsigned long);
#ifdef CORE_ENABLE_INLINES
  #include "Real.inl"
#else
  // friend functions for Real class
  // (need declarations in case they are not inlined)
  Real operator+(const Real&, const Real&);
  Real operator-(const Real&, const Real&);
  Real operator*(const Real&, const Real&);
  Real operator/(const Real&, const Real&);
  Real sqrt(const Real&);
  bool operator==(const Real&, const Real&);
  bool operator!=(const Real&, const Real&);
  bool operator< (const Real&, const Real&);
  bool operator<=(const Real&, const Real&);
  bool operator> (const Real&, const Real&);
  bool operator>=(const Real&, const Real&);
  int sign(const Real&);
  bool isZero(const Real&);
  int compare(const Real& e1, const Real& e2);
  BigInt floor(const Real&);
  BigInt ceil(const Real&);
  Real power(const Real&, unsigned long n);
  Real abs(const Real&);
  Real fabs(const Real&);
  std::ostream& operator<<(std::ostream&, const Real&);
  std::istream& operator>>(std::istream&, Real&);
#endif
  
#define CORE_REAL_ZERO Real::getZero()

CORE_END_NAMESPACE

#endif
