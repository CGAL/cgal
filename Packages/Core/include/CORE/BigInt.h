/******************************************************************
 * Core Library Version 1.6, June 2003
 * Copyright (c) 1995-2003 Exact Computation Project
 * 
 * File: BigInt.h
 *
 * Synopsis: a wrapper class of mpz in GMP
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

#ifndef CORE_BIGINT_H
#define CORE_BIGINT_H

#include <stdio.h> 
#include "gmp.h"
#include "CoreImpl.h"

CORE_BEGIN_NAMESPACE

/**************** Auxiliary classes ****************/

/* this is the same as gmp_allocated_string in gmp-impl.h
   since gmp-impl.h is not publicly available, I redefine it here
   I use a different name to avoid possible clashes */
struct __gmp_alloc_cstring
{
  char *str;
  __gmp_alloc_cstring(char *s) { str = s; }
  ~__gmp_alloc_cstring() { __gmp_free_func(str, strlen(str)+1); }
};

/// \class BigInt BigInt.h
/// \brief BigInt is a wrapper class of <tt>mpz</tt> in GMP
class BigInt
{
private:
  mpz_t mp;
public:
  /// \name Constructors and Destructor
  //@{
  /// default constructor
  BigInt() { mpz_init(mp); }
  /// copy constructor
  BigInt(const BigInt& z) { mpz_init_set(mp, z.mpz()); }

  /// constructor for <tt>signed char</tt>
  BigInt(signed char c) { mpz_init_set_si(mp, c); }
  /// constructor for <tt>unsigned char</tt>
  BigInt(unsigned char c) { mpz_init_set_ui(mp, c); } 

  /// constructor for <tt>signed int</tt>
  BigInt(signed int i) { mpz_init_set_si(mp, i); } 
  /// constructor for <tt>unsigned int</tt>
  BigInt(unsigned int i) { mpz_init_set_ui(mp, i); } 

  /// constructor for <tt>signed short int</tt>
  BigInt(signed short int s) { mpz_init_set_si(mp, s); } 
  /// constructor for <tt>unsigned short int </tt>
  BigInt(unsigned short int s) { mpz_init_set_ui(mp, s); } 

  /// constructor for <tt>signed long int</tt>
  BigInt(signed long int l) { mpz_init_set_si(mp, l); } 
  /// constructor for <tt>unsigned long int</tt>
  BigInt(unsigned long int l) { mpz_init_set_ui(mp, l); } 

  /// constructor for <tt>float</tt>
  BigInt(float f) { mpz_init_set_d(mp, f); } 
  /// constructor for <tt>double</tt>
  BigInt(double d) { mpz_init_set_d(mp, d); } 

  /// constructor for <tt>const char*</tt> with base
  BigInt(const char* s, int base = 0) { mpz_init_set_str(mp, s, base); } 
  /// constructor for <tt>std::string</tt> with base
  BigInt(const std::string& s, int base = 0) 
  { mpz_init_set_str(mp, s.c_str(), base); } 

  /// constructor for <tt>mpz_srcptr</tt>
  explicit BigInt(mpz_srcptr z) { mpz_init_set(mp, z); } 

  /// destructor
  ~BigInt() { mpz_clear(mp); }
  //@}

  /// \name Assignment Operators
  //@{ 
  /// = operator for <tt>BigInt</tt>
  BigInt& operator= (const BigInt& z) 
  { if (&z != this) mpz_set(mp, z.mpz()); return *this; } 
  
  /// = operator for <tt>signed char</tt> 
  BigInt& operator= (signed char c) { mpz_set_si(mp, c); return *this; } 
  /// = operator for <tt>unsigned char</tt> 
  BigInt& operator= (unsigned char c) { mpz_set_ui(mp, c); return *this; } 
  
  /// = operator for <tt>signed int</tt> 
  BigInt& operator= (signed int i) { mpz_set_si(mp, i); return *this; } 
  /// = operator for <tt>unsigned int</tt> 
  BigInt& operator= (unsigned int i) { mpz_set_ui(mp, i); return *this; } 
  
  /// = operator for <tt>signed short int</tt> 
  BigInt& operator= (signed short int s) { mpz_set_si(mp, s); return *this; } 
  /// = operator for <tt>unsigned short int</tt> 
  BigInt& operator= (unsigned short int s) { mpz_set_ui(mp, s); return *this; } 

  /// = operator for <tt>signed long int</tt> 
  BigInt& operator= (signed long int l) { mpz_set_si(mp, l); return *this; } 
  /// = operator for <tt>unsigned long int</tt> 
  BigInt& operator= (unsigned long int l) { mpz_set_ui(mp, l); return *this; } 
  
  /// = operator for <tt>float</tt>
  BigInt& operator= (float f) { mpz_set_d(mp, f); return *this; } 
  /// = operator for <tt>double</tt>
  BigInt& operator= (double d) { mpz_set_d(mp, d); return *this; } 
  
  /// = operator for <tt>const char*</tt> 
  BigInt& operator= (const char *s) { mpz_set_str(mp, s, 0); return *this; } 
  /// = operator for <tt>std::string</tt> 
  BigInt& operator= (const std::string &s) 
  { mpz_set_str(mp, s.c_str(), 0); return *this; }
  //@}

  /// \name Compound Assignment Operators
  //@{ 
  /// operator+= 
  BigInt& operator+= (const BigInt& z) 
  { mpz_add(mp, mp, z.mpz()); return *this; }
  /// operator-= 
  BigInt& operator-= (const BigInt& z) 
  { mpz_sub(mp, mp, z.mpz()); return *this; }
  /// operator*= 
  BigInt& operator*= (const BigInt& z) 
  { mpz_mul(mp, mp, z.mpz()); return *this; }
  /// operator/= 
  BigInt& operator/= (const BigInt& z) 
  { mpz_tdiv_q(mp, mp, z.mpz()); return *this; }
  /// operator%= 
  BigInt& operator%= (const BigInt& z) 
  { mpz_tdiv_r(mp, mp, z.mpz()); return *this; }
  /// operator%= 
  BigInt& operator&= (const BigInt& z) 
  { mpz_and(mp, mp, z.mpz()); return *this; }
  /// operator%= 
  BigInt& operator|= (const BigInt& z) 
  { mpz_ior(mp, mp, z.mpz()); return *this; }
  /// operator%= 
  BigInt& operator^= (const BigInt& z) 
  { mpz_xor(mp, mp, z.mpz()); return *this; }
  /// operator<<= 
  BigInt& operator<<= (long int l)
  { (l>=0)?mpz_mul_2exp(mp, mp, l):mpz_tdiv_q_2exp(mp, mp, -l); return *this; }
  /// operator>>=  
  BigInt& operator>>= (long int l)
  { (l>=0)?mpz_tdiv_q_2exp(mp, mp, l):mpz_mul_2exp(mp, mp, -l); return *this; }
  //@}
 
  /// \name Unary Minus, Increment, and Decrement Operators
  //@{
  /// unary minus
  BigInt operator- () const { BigInt c; mpz_neg(c.mpz(), mp); return c; }
  
  /// left increment operator (++i)
  BigInt& operator++ () { mpz_add_ui(mp, mp, 1); return *this; }
  /// left decrement operator (--i)
  BigInt& operator-- () { mpz_sub_ui(mp, mp, 1); return *this; }
  
  /// right increment operator (i++)
  BigInt operator++ (int) { BigInt a(*this); mpz_add_ui(mp, mp, 1); return a; }
  /// right deccrement operator (i--)
  BigInt operator-- (int) { BigInt a(*this); mpz_sub_ui(mp, mp, 1); return a; }
  //@}

  /// \name String Conversion Functions
  //@{
  /// set value from <tt>const char*</tt>
  /**
   * \returns 0 if the entire string is a valid number in base <tt>base</tt>. 
   *  Otherwise it returns -1
   */
  int fromString(const char* s, int base = 0)  
  { return mpz_set_str(mp, s, base); }
  /// convert to <tt>std::string</tt>
  std::string toString(int base = 10) const 
  { __gmp_alloc_cstring t(mpz_get_str(0, base, mp)); return std::string(t.str);}
  //@}

  /// \name Conversion Functions
  //@{
  /// return signed long int value
  signed long int longValue() const { return mpz_get_si(mp); }
  /// return unsigned long int value
  unsigned long int ulongValue() const { return mpz_get_ui(mp); }
  /// return double value
  double doubleValue() const { return mpz_get_d(mp); }   
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
  int sign() const { return mpz_sgn(mp); }
  /// absolute value function
  BigInt abs() const { BigInt c; mpz_abs(c.mpz(), mp); return c; }
  
  /// comparision function
  int cmp(const BigInt& q) const { return mpz_cmp(mp, q.mpz()); }
  /// absolute value comparision function
  int cmpabs(const BigInt& q) const { return mpz_cmpabs(mp, q.mpz()); }
  
  /// negation value function
  BigInt neg() const { BigInt c; mpz_neg(c.mpz(), mp); return c; }
  /// negate value function
  void negate() { mpz_neg(mp, mp); }
  
  /// bit length
  int bitLength() const { return mpz_sizeinbase(mp, 2); }
  
  /// get exponent of power 2
  unsigned long getBinExpo() const { return mpz_scan1(mp, 0); }
  /// get exponent of power k
  void getKaryExpo(BigInt& m, int& e, unsigned long k) const;

  /// get mpz pointer (const) (!!internal use!!)
  mpz_srcptr mpz() const { return mp; }
  /// get mpz pointer  (!!internal use!!)
  mpz_ptr mpz() { return mp; }
  //@}

  /// \name Type Checking functions
  //@{
  /// return true if it is even
  bool isEven() const { return mpz_even_p(mp); }
  /// return true if it is odd
  bool isOdd() const { return mpz_odd_p(mp); }
  /// return true if it fits a char
  bool isChar() const { return bitLength() <= 7; }
  /// return true if it fits an unsigned char 
  bool isUChar() const { return bitLength() <= 8 && sign() >= 0; }
  /// return true if it fits an int
  bool isInt() const { return mpz_fits_sint_p(mp) != 0; }
  /// return true if it fits an unsigned int 
  bool isUInt() const { return mpz_fits_uint_p(mp) != 0; }
  /// return true if it fits a short int
  bool isShort() const { return mpz_fits_sshort_p(mp) != 0; }
  /// return true if it fits an unsigned short int 
  bool isUShort() const { return mpz_fits_ushort_p(mp) != 0; }
  /// return true if it fits a long int
  bool isLong() const { return mpz_fits_slong_p(mp) != 0; }
  /// return true if it fits an unsigned long int 
  bool isULong() const { return mpz_fits_ulong_p(mp) != 0; }
  //@}
};

/// IO Stream operator<<
inline std::ostream& operator<<(std::ostream& o, const BigInt& a) 
{ return o << a.toString(); }
/// IO Stream operator>>
inline std::istream& operator>>(std::istream& i, BigInt& a) 
{ return ::operator>>(i, a.mpz()); }

/// operator+
inline BigInt operator+ (const BigInt& a, const BigInt& b) 
{ BigInt r(a); r += b; return r; }
/// operator-
inline BigInt operator- (const BigInt& a, const BigInt& b) 
{ BigInt r(a); r -= b; return r; }
/// operator*
inline BigInt operator* (const BigInt& a, const BigInt& b) 
{ BigInt r(a); r *= b; return r; }
/// operator/
inline BigInt operator/ (const BigInt& a, const BigInt& b) 
{ BigInt r(a); r /= b; return r; }
/// operator%
inline BigInt operator% (const BigInt& a, const BigInt& b) 
{ BigInt r(a); r %= b; return r; }
/// operator&
inline BigInt operator& (const BigInt& a, const BigInt& b) 
{ BigInt r(a); r &= b; return r; }
/// operator|
inline BigInt operator| (const BigInt& a, const BigInt& b) 
{ BigInt r(a); r |= b; return r; }
/// operator^
inline BigInt operator^ (const BigInt& a, const BigInt& b) 
{ BigInt r(a); r ^= b; return r; }
/// operator>>
inline BigInt operator>> (const BigInt& a, long int l)
{ BigInt r(a); r >>= l; return r; }
/// operator<<
inline BigInt operator<< (const BigInt& a, long int l)
{ BigInt r(a); r <<= l; return r; }

/// operator==
inline bool operator== (const BigInt& a, const BigInt& b)
{ return a.cmp(b) == 0; }
/// operator!=
inline bool operator!= (const BigInt& a, const BigInt& b) 
{ return a.cmp(b) != 0; }
/// operator>=
inline bool operator>= (const BigInt& a, const BigInt& b) 
{ return a.cmp(b) >= 0; }
/// operator>
inline bool operator> (const BigInt& a, const BigInt& b) 
{ return a.cmp(b) > 0; }
/// operator<=
inline bool operator<= (const BigInt& a, const BigInt& b) 
{ return a.cmp(b) <= 0; }
/// operator<
inline bool operator< (const BigInt& a, const BigInt& b) 
{ return a.cmp(b) < 0; }

/// sign
inline int sign(const BigInt& a) { return a.sign(); }
/// abs
inline BigInt abs(const BigInt& a) { return a.abs(); }
/// neg
inline BigInt neg(const BigInt& a) { return a.neg(); }
/// cmp
inline int cmp(const BigInt& a, const BigInt& b) { return a.cmp(b); }
/// cmpabs
inline int cmpabs(const BigInt& a, const BigInt& b) { return a.cmpabs(b); }

/// gcd
inline BigInt gcd(BigInt& a, BigInt& b) 
{ BigInt c; mpz_gcd(c.mpz(), a.mpz(), b.mpz()); return c; }
/// div_rem
inline void div_rem(BigInt& q, BigInt& r, const BigInt& a, const BigInt& b) 
{ mpz_tdiv_qr(q.mpz(), r.mpz(), a.mpz(), b.mpz()); }
/// power
inline void power(BigInt& c, const BigInt& a, unsigned long i) 
{ mpz_pow_ui(c.mpz(), a.mpz(), i); }
/// floorLg -- floor of log_2(a)
inline long floorLg(const BigInt& a) 
{ return ((a.sign() == 0) ? (-1) : (a.bitLength()-1)); }
/// ceilLg -- ceiling of log_2(a)
inline long ceilLg(const BigInt& a) {
  if (a.sign() == 0) return -1; 
  unsigned long len = a.bitLength();
  return (mpz_scan0(a.mpz(), 0) == len-1) ? (len-1) : len;
}

// return a gmp_randstate_t structure
extern gmp_randstate_t* getRandstate();
/// seed function
inline void seed(const BigInt& a)
{ BigInt tmp(a); gmp_randseed(*getRandstate(), tmp.mpz()); }
/// randomize function
inline BigInt randomize(const BigInt& a) 
{ BigInt c; mpz_urandomm(c.mpz(), *getRandstate(), a.mpz()); return c; }

CORE_END_NAMESPACE
#endif // CORE_BIGINT_H
