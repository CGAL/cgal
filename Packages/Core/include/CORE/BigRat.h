/******************************************************************
 * Core Library Version 1.6, June 2003
 * Copyright (c) 1995-2003 Exact Computation Project
 * 
 * File: BigRat.h
 *
 * Synopsis: a wrapper class of mpq in GMP
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


#ifndef CORE_BIGRAT_H
#define CORE_BIGRAT_H

#include "BigInt.h"

CORE_BEGIN_NAMESPACE

/// \class BigRat BigRat.h
/// \brief BigRat is a wrapper class of <tt>mpq</tt> in GMP
class BigRat
{
private:
  mpq_t mp;
public:
  /// \name Constructors and Destructor
  //@{
  /// default constructor
  BigRat() { mpq_init(mp); }
  /// copy constructor
  BigRat(const BigRat& q) { mpq_init(mp); mpq_set(mp, q.mpq()); }

  /// constructor for <tt>signed char</tt>
  BigRat(signed char c) { mpq_init(mp); mpq_set_si(mp, c, 1); }
  /// constructor for <tt>unsigned char</tt>
  BigRat(unsigned char c) { mpq_init(mp); mpq_set_ui(mp, c, 1); } 

  /// constructor for <tt>signed int</tt>
  BigRat(signed int i) { mpq_init(mp); mpq_set_si(mp, i, 1); } 
  /// constructor for <tt>unsigned int</tt>
  BigRat(unsigned int i) { mpq_init(mp); mpq_set_ui(mp, i, 1); } 

  /// constructor for <tt>signed short int</tt>
  BigRat(signed short int s) { mpq_init(mp); mpq_set_si(mp, s, 1); } 
  /// constructor for <tt>unsigned short int</tt>
  BigRat(unsigned short int s) { mpq_init(mp); mpq_set_ui(mp, s, 1); } 

  /// constructor for <tt>signed long int</tt>
  BigRat(signed long int l) { mpq_init(mp); mpq_set_si(mp, l, 1); } 
  /// constructor for <tt>unsigned long int</tt>
  BigRat(unsigned long int l) { mpq_init(mp); mpq_set_ui(mp, l, 1); } 

  /// constructor for <tt>float</tt>
  BigRat(float f) { mpq_init(mp); mpq_set_d(mp, f); } 
  /// constructor for <tt>double</tt>
  BigRat(double d) { mpq_init(mp); mpq_set_d(mp, d); } 

  /// constructor for <tt>const char*</tt> with base
  BigRat(const char *s, int base = 0)
  { mpq_init(mp); mpq_set_str(mp, s, base); } 
  /// constructor for <tt>std::string</tt> with base
  BigRat(const std::string &s, int base = 0) 
  { mpq_init(mp); mpq_set_str(mp, s.c_str(), base); } 

  /// constructor for <tt>mpq_srcptr</tt>
  explicit BigRat(mpq_srcptr q) { mpq_init(mp); mpq_set(mp, q); } 

  /// constructor for <tt>BigInt</tt>
  BigRat(const BigInt& z) 
  { mpq_init(mp); mpq_set_z(mp, z.mpz()); }
  /// constructor for two <tt>BigInts</tt> 
  BigRat(const BigInt& n, const BigInt& d) { 
    mpq_init(mp); 
    mpz_set(mpq_numref(mp), n.mpz()); 
    mpz_set(mpq_denref(mp), d.mpz());
    mpq_canonicalize(mp);
  }

  /// destructor
  ~BigRat() { mpq_clear(mp); }
  //@}

  /// \name Assignment Operators
  //@{ 
  /// = operator for <tt>BigRat</tt>
  BigRat& operator= (const BigRat& q) 
  { if (&q != this) mpq_set(mp, q.mpq()); return *this; } 
  
  /// = operator for <tt>signed char</tt> 
  BigRat& operator= (signed char c) { mpq_set_si(mp, c, 1); return *this; } 
  /// = operator for <tt>unsigned char</tt> 
  BigRat& operator= (unsigned char c) { mpq_set_ui(mp, c, 1); return *this; } 
  
  /// = operator for <tt>signed int</tt> 
  BigRat& operator= (signed int i) { mpq_set_si(mp, i, 1); return *this; } 
  /// = operator for <tt>unsigned int</tt> 
  BigRat& operator= (unsigned int i) { mpq_set_ui(mp, i, 1); return *this; } 
  
  /// = operator for <tt>signed short int</tt> 
  BigRat& operator= (signed short int s) { mpq_set_si(mp, s, 1); return *this; } 
  /// = operator for <tt>unsigned short int</tt> 
  BigRat& operator= (unsigned short int s) 
  { mpq_set_ui(mp, s, 1); return *this; } 
  
  /// = operator for <tt>signed long int</tt> 
  BigRat& operator= (signed long int l) { mpq_set_si(mp, l, 1); return *this; } 
  /// = operator for <tt>unsigned long int</tt> 
  BigRat& operator= (unsigned long int l) 
  { mpq_set_ui(mp, l, 1); return *this; } 
  
  /// = operator for <tt>float</tt>
  BigRat& operator= (float f) { mpq_set_d(mp, f); return *this; } 
  /// = operator for <tt>double</tt>
  BigRat& operator= (double d) { mpq_set_d(mp, d); return *this; } 
  
  /// = operator for <tt>const char*</tt> 
  BigRat& operator= (const char *s) { mpq_set_str(mp, s, 0); return *this; } 
  /// = operator for <tt>std::string</tt> 
  BigRat& operator= (const std::string &s) 
  { mpq_set_str(mp, s.c_str(), 0); return *this; }
  //@}

  /// \name Compound Assignment Operators
  //@{ 
  /// operator+= 
  BigRat& operator+= (const BigRat& q) 
  { mpq_add(mp, mp, q.mpq()); return *this; }
  /// operator-= 
  BigRat& operator-= (const BigRat& q) 
  { mpq_sub(mp, mp, q.mpq()); return *this; }
  /// operator*= 
  BigRat& operator*= (const BigRat& q) 
  { mpq_mul(mp, mp, q.mpq()); return *this; }
  /// operator/= 
  BigRat& operator/= (const BigRat& q) 
  { mpq_div(mp, mp, q.mpq()); return *this; }
  /// operator<<= 
  BigRat& operator<<= (long int l)
  { (l>=0)?mpq_mul_2exp(mp, mp, l):mpq_div_2exp(mp, mp, -l); return *this; }
  /// operator>>=  
  BigRat& operator>>= (long int l)
  { (l>=0)?mpq_div_2exp(mp, mp, l):mpq_mul_2exp(mp, mp, -l); return *this; }
  //@}
 
  /// \name Unary Minus, Increment, and Decrement Operators
  //@{
  /// unary minus
  BigRat operator- () const { BigRat c; mpq_neg(c.mpq(), mp); return c; }
  
  /// left increment operator (++i)
  BigRat& operator++ () 
  { mpz_add(mpq_numref(mp), mpq_numref(mp), mpq_denref(mp)); return *this; }
  /// left decrement operator (--i)
  BigRat& operator-- () 
  { mpz_sub(mpq_numref(mp), mpq_numref(mp), mpq_denref(mp)); return *this; }
  
  /// right increment operator (i++)
  BigRat operator++ (int) {
    BigRat a(*this);
    mpz_add(mpq_numref(mp), mpq_numref(mp), mpq_denref(mp));
    return a;
  }
  /// right deccrement operator (i--)
  BigRat operator-- (int) {
    BigRat a(*this);
    mpz_add(mpq_numref(mp), mpq_numref(mp), mpq_denref(mp));
    return a;
  }
  //@}

  /// \name String Conversion Functions
  //@{
  /// set value from <tt>const char*</tt>
  /**
   * \returns 0 if the entire string is a valid number in base <tt>base</tt>. 
   *  Otherwise it returns -1
   */
  int fromString(const char* s, int base = 0)  
  { return mpq_set_str(mp, s, base); }
  /// convert to <tt>std::string</tt>
  std::string toString(int base = 10) const 
  { __gmp_alloc_cstring t(mpq_get_str(0, base, mp)); return std::string(t.str);}
  //@}

  /// \name Conversion Functions
  //@{
  /// return double value
  double doubleValue() const { return mpq_get_d(mp); }   
  /// return BigInt value
  BigInt BigIntValue() const { return numerator()/denominator(); }
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
  /// Canonicalize 
  void canonicalize() { mpq_canonicalize(mp); }

  /// return sign
  int sign() const { return mpq_sgn(mp); }
  /// return absolute value
  BigRat abs() const { BigRat c; mpq_abs(c.mpq(), mp); return c; }
  
  /// return negation value
  BigRat neg() const { BigRat c; mpq_neg(c.mpq(), mp); return c; }
  /// negate value
  void negate() { mpq_neg(mp, mp); }
  
  /// comparision function
  int cmp(const BigRat& q) const { return mpq_cmp(mp, q.mpq()); }
  
  /// return the numerater (const)
  const BigInt numerator() const { return BigInt(mpq_numref(mp)); }
  /// return the numerater
  BigInt numerator() { return BigInt(mpq_numref(mp)); }
  
  /// return the denominator (const)
  const BigInt denominator() const { return BigInt(mpq_denref(mp)); }
  /// return the denominator
  BigInt denominator() { return BigInt(mpq_denref(mp)); }
  
  /// return mpq pointer (const)  (!!internal use!!)
  mpq_srcptr mpq() const { return mp; }
  /// return mpq pointer  (!!internal use!!)
  mpq_ptr mpq() { return mp; }  
  
  /// return mpz pointer of numerator (const)  (!!internal use!!)
  mpz_srcptr num_mpz() const { return mpq_numref(mp); }
  /// return mpz pointer of numerator  (!!internal use!!)
  mpz_ptr num_mpz() { return mpq_numref(mp); }
  /// return mpz pointer of denominator (const)  (!!internal use!!)
  
  mpz_srcptr den_mpz() const { return mpq_denref(mp); }
  /// return mpz pointer of denominator  (!!internal use!!)
  mpz_ptr den_mpz() { return mpq_denref(mp); }
  //@}
};

/// IO stream operator <<
inline std::ostream& operator<<(std::ostream& o, const BigRat& a) 
{ return o << a.toString(); }
/// IO stream operator >>
inline std::istream& operator>>(std::istream& i, BigRat& a) 
{ return ::operator>>(i, a.mpq()); }

/// operator+
inline BigRat operator+ (const BigRat& a, const BigRat& b) 
{ BigRat r(a); r += b; return r; }
/// operator-
inline BigRat operator- (const BigRat& a, const BigRat& b) 
{ BigRat r(a); r -= b; return r; }
/// operator*
inline BigRat operator* (const BigRat& a, const BigRat& b) 
{ BigRat r(a); r *= b; return r; }
/// operator/
inline BigRat operator/ (const BigRat& a, const BigRat& b) 
{ BigRat r(a); r /= b; return r; }
/// operator>>
inline BigRat operator>> (const BigRat& a, long int l)
{ BigRat r(a); r >>= l; return r; }
/// operator<<
inline BigRat operator<< (const BigRat& a, long int l)
{ BigRat r(a); r <<= l; return r; }

/// operator==
inline bool operator== (const BigRat& a, const BigRat& b) 
{ return a.cmp(b) == 0; }
/// operator!=
inline bool operator!= (const BigRat& a, const BigRat& b) 
{ return a.cmp(b) != 0; }
/// operator>=
inline bool operator>= (const BigRat& a, const BigRat& b) 
{ return a.cmp(b) >= 0; }
/// operator>
inline bool operator> (const BigRat& a, const BigRat& b) 
{ return a.cmp(b) > 0; }
/// operator<=
inline bool operator<= (const BigRat& a, const BigRat& b) 
{ return a.cmp(b) <= 0; }
/// operator<
inline bool operator< (const BigRat& a, const BigRat& b) 
{ return a.cmp(b) < 0; }

/// sign
inline int sign(const BigRat& a) { return a.sign(); }
/// abs
inline BigRat abs(const BigRat& a) { return a.abs(); }
/// neg
inline BigRat neg(const BigRat& a) { return a.neg(); }
/// cmp
inline int cmp(const BigRat& a, const BigRat& b) { return a.cmp(b); }

CORE_END_NAMESPACE
#endif // CORE_BIGRAT_H

