/******************************************************************
 * Core Library Version 1.5, August 2002
 * Copyright (c) 1995-2002 Exact Computation Project
 * 
 * File: BigInt.h
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

#include "CoreImpl.h"

/*****************************************************************
 * The Following is adapted from LiDIA's bigint.h and bigint_def.h
 *****************************************************************/

#include <gmp.h>

CORE_BEGIN_NAMESPACE

//typedef mp_limb_t UWtype;

#ifndef BITS_PER_CHAR
#define BITS_PER_CHAR	8
#endif

// class bigint: modified from LiDia, Zilin Du, May 30, 2001

class bigint
{

  /**
  ** the C type we use to represent a bigint
  **/
public:
  mpz_t I;
private:

  /**
   **  input / output  utilities, implemented in bigint_share.c
   **/

  static int chars_per_line;    // number of characters in one line
                                // when printing a bigint


  // The next 5 functions are implemented in bigint_share.c
  // and exist only once for all interfaces.

  static void allocate (char * &s, int old_size, int new_size);

public:
  static void append_char (char * &s, int& sz, int pos, char c);
  static int  skip_backslash_new_line (std::istream & in);

  void scan  (std::istream & in);
  void print (std::ostream & out, char *s) const;
  

 public:

  // The next 2 functions are implemented in bigint_share.c
  // and exist only once for all interfaces.

  static void set_chars_per_line (int cpl);
  static int  get_chars_per_line ();


  /**
  ** constructors and destructor; we could leave out some of these
  **/

  bigint();
  bigint(int i);
  bigint(unsigned int ui);
  bigint(long l);
  bigint(unsigned long ul);
  bigint(const char* str); // new constructor
  bigint(const mpz_t & a);
  bigint(const bigint & a);
  ~bigint();

  /**
  ** inline member functions
  **/
  
  int bit(unsigned int i) const;
  int length() const;
  int bit_length() const;
  int sign() const;

  bool is_odd() const;
  bool is_even() const;

  friend bool is_odd  (const bigint & a);
  friend bool is_even (const bigint & a);

#ifndef HEADBANGER

  bool is_positive() const;
  bool is_negative() const;
  bool is_zero() const;
  bool is_gt_zero() const;
  bool is_ge_zero() const;
  bool is_lt_zero() const;
  bool is_le_zero() const;
  bool is_one() const;

  friend bool is_positive (const bigint & a);
  friend bool is_negative (const bigint & a);
  friend bool is_zero (const bigint & a);
  friend bool is_one (const bigint & a);

  bool intify(int &i) const;
  bool longify(long &i) const;
  int abs_compare(const bigint & a) const;
  int compare(const bigint & a) const;
  unsigned long most_significant_digit() const;
  unsigned long least_significant_digit() const;

  bigint characteristic () const;

  /**
  ** the next two definitions are needed by bigfloat
  **/

  static const double radix();
  static const int bits_per_digit();

  void absolute_value();
  void abs();
  void negate();
  void assign_zero();
  void assign_one();
  void assign(int i);
  void assign(long ui);
  void assign(unsigned long ui);
  void assign(double);
  void assign(const bigint & a);
  void multiply_by_2();
  void divide_by_2();

  /**
  ** type checking
  **/

  friend bool is_char(const bigint & a);
  friend bool is_uchar(const bigint & a);
  friend bool is_short(const bigint & a);
  friend bool is_ushort(const bigint & a);
  friend bool is_int(const bigint & a);
  friend bool is_uint(const bigint & a);
  friend bool is_long(const bigint & a);
  friend bool is_ulong(const bigint & a);

#endif

  /**
  ** assignments
  **/

  int operator = (int i);
  long operator = (long l);
  unsigned long operator = (unsigned long ul);
  double operator = (double d);
  bigint & operator = (const bigint & a);

  /**
  ** comparisons
  **/

  friend bool operator == (const bigint & a, const bigint & b);
  friend bool operator != (const bigint & a, const bigint & b);
  friend bool operator > (const bigint & a, const bigint & b);
  friend bool operator >= (const bigint & a, const bigint & b);
  friend bool operator < (const bigint & a, const bigint & b);
  friend bool operator <= (const bigint & a, const bigint & b);

  /**
  ** operator overloading
  **/

  friend bigint operator - (const bigint & a);
  friend bigint operator + (const bigint & a, const bigint & b);
  friend bigint operator - (const bigint & a, const bigint & b);
  friend bigint operator *(const bigint & a, const bigint & b);
  friend bigint operator / (const bigint & a, const bigint & b);
  friend bigint operator % (const bigint & a, const bigint & b);
  friend bigint operator << (const bigint & a, long u);
  friend bigint operator >> (const bigint & a, long u);
  friend bigint operator & (const bigint & a, const bigint & b);
  friend bigint operator | (const bigint & a, const bigint & b);
  friend bigint operator ^ (const bigint & a, const bigint & b);

  bigint & operator += (const bigint & a);
  bigint & operator -= (const bigint & a);
  bigint & operator *= (const bigint & a);
  bigint & operator /= (const bigint & a);
  bigint & operator %= (const bigint & a);
  bigint & operator <<= (long ui);
  bigint & operator >>= (long ui);
  bigint & operator &= (const bigint & a);
  bigint & operator |= (const bigint & a);
  bigint & operator ^= (const bigint & a);

  bigint operator~ () const;
  bigint & operator++ ();
  bigint & operator-- ();
  bigint   operator++ (int);
  bigint   operator-- (int);
  int operator ! () const;
  
  unsigned long getBinExpo() const { return mpz_scan1(I, 0); }
  void getKaryExpo(bigint& m, int& e, unsigned long k) const;

#ifndef HEADBANGER

  /**
  ** Procedural versions
  **/

  friend void lidia_negate(bigint & a, const bigint & b);
  friend void add(bigint & c, const bigint & a, const bigint & b);
  friend void subtract(bigint & c, const bigint & a, const bigint & b);
  friend void multiply(bigint & c, const bigint & a, const bigint & b);
  friend void divide(bigint & c, const bigint & a, const bigint & b);
  friend void remainder(bigint & c, const bigint & a, const bigint & b);
  friend void div_rem(bigint & q, bigint & r, const bigint & a, const bigint & b);
  friend void shift_left(bigint & c, const bigint & a, long ui);
  friend void shift_right(bigint & c, const bigint & a, long ui);
  friend void power(bigint & c, const bigint & a, const bigint & b);
  friend void power(bigint & c, const bigint & a, long i);
  //friend void and(bigint & c, const bigint & a, const bigint & b);
  //friend void or(bigint & c, const bigint & a, const bigint & b);
  //friend void xor(bigint & c, const bigint & a, const bigint & b);
  //friend void not(bigint & c, const bigint & a);
  friend void inc(bigint & c);
  friend void dec(bigint & c);

  friend void add(bigint & c, const bigint & a, long i);
  friend void subtract(bigint & c, const bigint & a, long i);
  friend void multiply(bigint & c, const bigint & a, long i);
  friend void divide(bigint & c, const bigint & a, long i);
  friend void remainder(long &r, const bigint & a, long i);
  friend void remainder(unsigned long &r, const bigint & a, unsigned long i);
  friend long remainder(const bigint & a, long i);
  friend void div_rem(bigint & q, long &r, const bigint & a, long i);
  friend void invert(bigint & a, const bigint & b);

#endif

  /**
  ** gcd's
  **/

  friend bigint gcd(const bigint & a, const bigint & b);
  friend bigint bgcd(const bigint & a, const bigint & b);
  friend bigint dgcd(const bigint & a, const bigint & b);
  friend bigint xgcd(bigint & u, bigint & v, const bigint & a, const bigint & b);
  friend bigint xgcd_left(bigint & u, const bigint & a, const bigint & b);
  friend bigint xgcd_right(bigint & v, const bigint & a, const bigint & b);

  /**
  ** functions
  **/

  friend bigint abs(const bigint & a);
  friend void seed(const bigint & a);
  friend bigint randomize(const bigint & a);
  void randomize(const bigint & a);
  friend double dbl(const bigint & a);
//  friend xdouble xdbl(const bigint & a);

  friend void sqrt(bigint & a, const bigint & b);
  friend void square(bigint & a, const bigint & b);
  friend void swap(bigint & a, bigint & b);

  /**
  ** input / output
  **/

  friend std::istream & operator >> (std::istream & in, bigint & a);
  friend std::ostream & operator << (std::ostream & out, const bigint & a);

  friend int string_to_bigint(const char *s, bigint & a, int base = 10);
  friend int bigint_to_string(const bigint & a, char *s, int base = 10);

  // Note: Remove below 4 functions since it seems we never use it in Core Library
  //    Zilin Du, June 14, 2001
  
  /**
  ** using fread/fwrite
  **/
  
  //void read_from_file(FILE * fp);
  //void write_to_file(FILE * fp);

  /**
  ** using fscanf/fprintf
  **/
  
  //void scan_from_file(FILE * fp);
  //void print_to_file(FILE * fp);
  
  // Note: Add 2 function to implement read/write BigInt in our own format
  //       Zilin Du, June 14, 2001
  void read_from_file(std::istream& in, long maxLength = 0);
  void write_to_file(std::ostream& out, int base = 10, int charsPerLine = 80);
  static int skip_comment_line (std::istream & in);
  static void read_string(std::istream& in, char* &buffer, int sz);
  static void read_base_number(std::istream& in, bigint& m, long bits, long maxBits);
  static void write_base_number(std::ostream& out, char* buffer, 
                        int length, int base, int charsPerLine);
  
  static const bigint& getAnonymous();
};

#ifdef CORE_ENABLE_INLINES
  #include "BigInt.inl"
#else
CORE_INLINE void lidia_error_handler(const char *f, const char *m);
  bool is_odd  (const bigint & a);
  bool is_even (const bigint & a);
  bool is_positive (const bigint & a);
  bool is_negative (const bigint & a);
  bool is_zero (const bigint & a);
  bool is_one (const bigint & a);
  bool is_char(const bigint & a);
  bool is_uchar(const bigint & a);
  bool is_short(const bigint & a);
  bool is_ushort(const bigint & a);
  bool is_int(const bigint & a);
  bool is_uint(const bigint & a);
  bool is_long(const bigint & a);
  bool is_ulong(const bigint & a);
  bool operator == (const bigint & a, const bigint & b);
  bool operator != (const bigint & a, const bigint & b);
  bool operator > (const bigint & a, const bigint & b);
  bool operator >= (const bigint & a, const bigint & b);
  bool operator < (const bigint & a, const bigint & b);
  bool operator <= (const bigint & a, const bigint & b);
  bigint operator - (const bigint & a);
  bigint operator + (const bigint & a, const bigint & b);
  bigint operator - (const bigint & a, const bigint & b);
  bigint operator *(const bigint & a, const bigint & b);
  bigint operator / (const bigint & a, const bigint & b);
  bigint operator % (const bigint & a, const bigint & b);
  bigint operator << (const bigint & a, long u);
  bigint operator >> (const bigint & a, long u);
  bigint operator & (const bigint & a, const bigint & b);
  bigint operator | (const bigint & a, const bigint & b);
  bigint operator ^ (const bigint & a, const bigint & b);
  void lidia_negate(bigint & a, const bigint & b);
  void add(bigint & c, const bigint & a, const bigint & b);
  void subtract(bigint & c, const bigint & a, const bigint & b);
  void multiply(bigint & c, const bigint & a, const bigint & b);
  void divide(bigint & c, const bigint & a, const bigint & b);
  void remainder(bigint & c, const bigint & a, const bigint & b);
  void div_rem(bigint & q, bigint & r, const bigint & a, const bigint & b);
  void shift_left(bigint & c, const bigint & a, long ui);
  void shift_right(bigint & c, const bigint & a, long ui);
  void power(bigint & c, const bigint & a, const bigint & b);
  void power(bigint & c, const bigint & a, long i);
  void and(bigint & c, const bigint & a, const bigint & b);
  void or(bigint & c, const bigint & a, const bigint & b);
  void xor(bigint & c, const bigint & a, const bigint & b);
  void not(bigint & c, const bigint & a);
  void inc(bigint & c);
  void dec(bigint & c);

  void add(bigint & c, const bigint & a, long i);
  void subtract(bigint & c, const bigint & a, long i);
  void multiply(bigint & c, const bigint & a, long i);
  void divide(bigint & c, const bigint & a, long i);
  void remainder(long &r, const bigint & a, long i);
  void remainder(unsigned long &r, const bigint & a, unsigned long i);
  long remainder(const bigint & a, long i);
  void div_rem(bigint & q, long &r, const bigint & a, long i);
  void invert(bigint & a, const bigint & b);
  bigint gcd(const bigint & a, const bigint & b);
  bigint bgcd(const bigint & a, const bigint & b);
  bigint dgcd(const bigint & a, const bigint & b);
  bigint xgcd(bigint & u, bigint & v, const bigint & a, const bigint & b);
  bigint xgcd_left(bigint & u, const bigint & a, const bigint & b);
  bigint xgcd_right(bigint & v, const bigint & a, const bigint & b);
  bigint abs(const bigint & a);
  void seed(const bigint & a);
  bigint randomize(const bigint & a);
  double dbl(const bigint & a);

  void sqrt(bigint & a, const bigint & b);
  void square(bigint & a, const bigint & b);
  void swap(bigint & a, bigint & b);

  std::istream & operator >> (std::istream & in, bigint & a);
  std::ostream & operator << (std::ostream & out, const bigint & a);

  int string_to_bigint(const char *s, bigint & a, int base);
  int bigint_to_string(const bigint & a, char *s, int base);

  long bigIntToLong(const bigint & a);
  double bigIntToDouble(const bigint & a);
  int lg(const bigint & a);
  int ceilLg(const bigint & a);
  int floorLg(const bigint & a);
  int compare(const bigint &a, const bigint &b);
  int sign(const bigint &a);  

#endif

typedef bigint BigInt;

#define CORE_BIGINT_ZERO bigint::getAnonymous()

CORE_END_NAMESPACE

#endif
