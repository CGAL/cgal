/******************************************************************
 * Core Library Version 1.5, August 2002
 * Copyright (c) 1995-2002 Exact Computation Project
 * 
 * File: BigRat.h
 * Synopsis:
 *      This is Core Library's big rational class. 
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

/*****************************************************************
 * Following is from LiDIA's bigrational.h.   Zilin Du, May 30, 2001
 *****************************************************************/
#include "CoreImpl.h"
#include "BigInt.h"

CORE_BEGIN_NAMESPACE

void newton_root(bigint & b, const bigint & a, int n);

class bigrational
{
  /**
   ** the C++ type we use to represent a bigrational
   **/

  bigint num, den;
  void normalize();

 public:

  /**
   ** constructors and destructor we could leave out some of these
   **/

  bigrational();
  bigrational(int i);
  bigrational(long l);
  bigrational(unsigned long ul);
  bigrational(double d);
  bigrational(const char* s);
  bigrational(const bigint & n);
  bigrational(const bigint & n, const bigint & d);
  bigrational(const bigrational & a);
  ~bigrational();

#ifndef HEADBANGER

  /**
   ** inline member functions
   **/

  int sign() const;
  bool is_positive() const;
  bool is_negative() const;
  bool is_zero() const;
  bool is_gt_zero() const;
  bool is_ge_zero() const;
  bool is_lt_zero() const;
  bool is_le_zero() const;
  bool is_one() const;
  bool intify(int &i) const;
  bool longify(long &i) const;
  int abs_compare(const bigrational & a) const;
  int compare(const bigrational & a) const;

  void absolute_value();
  void negate();
  void invert();

  bigint numerator() const;
  bigint denominator() const;

  void multiply_by_denominator();

  void assign_zero();
  void assign_one();
  void assign(const bigint & n, const bigint & d);
  void assign(const bigrational & a);
  void multiply_by_2();
  void divide_by_2();

  // for use in elliptic curve code
  bigint characteristic() const;

#endif

  /**
   ** type checking
   **/

friend bool is_bigint(const bigrational & a);

  /**
   ** assignments
   **/

  int operator = (int i);
  long operator = (long l);
  unsigned long operator = (unsigned long ul);
  double operator = (double d);
  bigint operator = (const bigint & a);
  bigrational & operator = (const bigrational & a);

  /**
   ** comparisons
   **/

friend bool operator == (const bigrational & a, const bigrational & b);
friend bool operator != (const bigrational & a, const bigrational & b);
friend bool operator > (const bigrational & a, const bigrational & b);
friend bool operator >= (const bigrational & a, const bigrational & b);
friend bool operator < (const bigrational & a, const bigrational & b);
friend bool operator <= (const bigrational & a, const bigrational & b);

  /**
   ** operator overloading
   **/

friend bigrational operator - (const bigrational & a);
friend bigrational operator + (const bigrational & a, const bigrational & b);
friend bigrational operator - (const bigrational & a, const bigrational & b);
friend bigrational operator *(const bigrational & a, const bigrational & b);
friend bigrational operator / (const bigrational & a, const bigrational & b);
friend bigrational operator << (const bigrational & a, unsigned long ui);
friend bigrational operator >> (const bigrational & a, unsigned long ui);

  bigrational & operator += (const bigrational & a);
  bigrational & operator -= (const bigrational & a);
  bigrational & operator *= (const bigrational & a);
  bigrational & operator /= (const bigrational & a);
  bigrational & operator <<= (unsigned long ui);
  bigrational & operator >>= (unsigned long ui);

  bigrational & operator++ ();
  bigrational & operator-- ();
  int operator ! ();

#ifndef HEADBANGER

  /**
   ** Procedural versions
   **/

friend void invert(bigrational & c, const bigrational & a);
friend void lidia_negate(bigrational & a, const bigrational & b);
friend void add(bigrational & c, const bigrational & a, const bigrational & b);
friend void subtract(bigrational & c, const bigrational & a, const bigrational & b);
friend void multiply(bigrational & c, const bigrational & a, const bigrational & b);
friend void divide(bigrational & c, const bigrational & a, const bigrational & b);
friend void shift_left(bigrational & c, const bigrational & a, long ui);
friend void shift_right(bigrational & c, const bigrational & a,long ui);
friend void power(bigrational & c, const bigrational & a, const bigint & b);
friend void power(bigrational & c, const bigrational & a, long i);
friend void inc(bigrational & c);
friend void dec(bigrational & c);

friend void add(bigrational & c, const bigrational & a, long i);
friend void subtract(bigrational & c, const bigrational & a, long i);
friend void multiply(bigrational & c, const bigrational & a, long i);
friend void divide(bigrational & c, const bigrational & a, long i);

friend void add(bigrational & c, long i, const bigrational & a);
friend void subtract(bigrational & c, long i, const bigrational & a);
friend void multiply(bigrational & c, long i, const bigrational & a);
friend void divide(bigrational & c, long i, const bigrational & a);

#endif

  /**
   ** functions
   **/

friend bigrational abs(const bigrational & a);
friend bigrational inverse(const bigrational & a);

friend bigint numerator(const bigrational & a);
friend bigint denominator(const bigrational & a);
friend bigint round(const bigrational & a);
friend bigint floor(const bigrational & a);
friend bigint ceiling(const bigrational & a);
friend bigint truncate(const bigrational & a);

friend double dbl(const bigrational & a);

friend void square(bigrational & a, const bigrational & b);
friend void swap(bigrational & a, bigrational & b);

  /**
   ** input / output
   **/

friend std::istream & operator >> (std::istream & in, bigrational & a);
friend std::ostream & operator << (std::ostream & out, const bigrational & a);

friend int string_to_bigrational(char *s, char *t, bigrational & a);
friend int bigrational_to_string(const bigrational & a, char *s, char *t);

// Note: Remove below 4 functions since it seems we never use it in Core Library
//      Zilin Du, June 14, 2001
   /**
   ** using fread/fwrite
   **/

//  void read_from_file(FILE * fp);
//  void write_to_file(FILE * fp);

  /**
   ** using fscanf/fprintf
   **/

//  void scan_from_file(FILE * fp);
//  void print_to_file(FILE * fp);

};

bool square_root(bigrational & root, const bigrational& s);
bool cube_root(bigrational & root, const bigrational& s);

#ifdef CORE_ENABLE_INLINES
#include "BigRat.inl"
#endif

typedef bigrational BigRat;

CORE_END_NAMESPACE
#endif

