/******************************************************************
 * Core Library Version 1.5, August 2002
 * Copyright (c) 1995-2002 Exact Computation Project
 * 
 * File: BigRat.cpp
 * 
 * Synopsis:
 *      Provides the big rational numbers in Core Library
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
//
// LiDIA - a library for computational number theory
//   Copyright (c) 1994, 1995 by the LiDIA Group
//
// File        : bigrational.c 
// Author      : Thomas Papanikolaou (TP)
// Last change : TP, Jan 29 1995, initial version 
//               VM, Feb 7 1995, modified operators +, -, *, /
//                   to use the procedures add(), subtract(), ...
//                   Added a space ' ' in the output routine.
//               MM, Sep 11 1996, added 
//                   assign(const bigint&,const bigint&)
//               TPf, Sep 24 1999, corrected bug in ceiling (reported by 
//                   Georg Kux)
//               MM, Oct 18 1999, added
//                   characteristic() function

#include "BigRat.h"

CORE_BEGIN_NAMESPACE

#ifndef CORE_ENABLE_INLINES
#include "BigRat.inl"
#endif

typedef mp_limb_t base_digit;

bigrational::bigrational(double d)
{
  int exponent, sign = (d < 0.0);
  double fraction;
  base_digit m;
  long ndigits = 0;

  num.assign_zero();
  den.assign_one();

  if (d)
  {
    if (sign)
      d = -d;

    fraction = frexp(d, &exponent);

    while (fraction != 0)
    {
      shift_left(num, num, bigint::bits_per_digit()); 
      fraction = ldexp(fraction, bigint::bits_per_digit());
      m = (base_digit) fraction;
      fraction -= m;
      ndigits++;
      add(num, num, bigint(m));
    }
    // the product on the right side is never > MAXINT
    exponent -= (int)(ndigits * bigint::bits_per_digit());
    if (exponent < 0)
      shift_left(den, den, -exponent);
    else
      shift_left(num, num, exponent);
    if (sign)
      num.negate();
  }
  this->normalize();
}

bigrational::bigrational(const char* str)
{
  int len = strlen(str);
  int i=0; 
  const char* p = str;
  while ( (i<len) && (*p) != '/' ) {
    p++; i++;
  }
  if (i < len) {
    char* numStr = new char[i+1];
    strncpy(numStr, str, i); numStr[i]='\0';
    BigInt n(numStr);
    //std::cout << "num=|" << numStr << "|" << std::endl;
    char* denStr = new char[len-i];
    strncpy(denStr, (p+1), len-i-1); denStr[len-i-1]='\0';
    //std::cout << "den=|" << denStr << "|" << std::endl;
    BigInt d(denStr);
    num.assign(n);
    den.assign(d);
    this->normalize();
    //std::cout << "num=" << num << ", den=" << den << std::endl;
  } else {
    BigInt m(str);
    num.assign(m);
    den.assign_one();
  }
}

double bigrational::operator = (double d)
{
  int exponent, sign = (d < 0.0);
  double fraction;
  base_digit m;
  long ndigits = 0;

  num.assign_zero();
  den.assign_one();

  if (d)
  {
    if (sign)
      d = -d;

    fraction = frexp(d, &exponent);

    while (fraction != 0)
    {
      shift_left(num, num, bigint::bits_per_digit()); 
      fraction = ldexp(fraction, bigint::bits_per_digit());
      m = (base_digit) fraction;
      fraction -= m;
      ndigits++;
      add(num, num, bigint(m));
    }
    // the product on the right side is never > MAXINT
    exponent -= (int)(ndigits * bigint::bits_per_digit());
    if (exponent < 0)
      shift_left(den, den, -exponent);
    else
      shift_left(num, num, exponent);
    if (sign)
      num.negate();
  }
  this->normalize();
  return d;
}

void
add(bigrational & c, const bigrational & a, const bigrational & b)
{
  if (a.is_zero())
  {
    c.num.assign(b.num);
    c.den.assign(b.den);
  }
  else if (b.is_zero())
  {
    c.num.assign(a.num);
    c.den.assign(a.den);
  }
  else
  {
    bigint g = gcd(a.den, b.den), h;
    if (g.is_one())
    {
      // c.num   a.num * b.den + a.den * b.num;
      // ----- = -----------------------------
      // c.den         a.den * b.den;

      multiply(h, a.num, b.den);
      multiply(c.num, a.den, b.num);
      add(c.num, c.num, h);
      multiply(c.den, a.den, b.den);
    }
    else
    {
      // bigint t = a.num * (b.den / g) + b.num * (a.den / g);
      // bigint h = gcd(t, g);
      // c.num = t / h;
      // c.den = (a.den / g) * (b.den / h);

      bigint t, s, ss;
      divide(s, b.den, g);
      multiply(t, a.num, s);
      divide(s, a.den, g);
      multiply(ss, s, b.num);
      add(t, t, ss);
      h = gcd(t, g);
      divide(c.num, t, h);
      divide(t, b.den, h);
      multiply(c.den, s, t);
    }
    if (c.den.is_negative())
    {
      c.num.negate();
      c.den.negate();
    }
  }
}

void
subtract(bigrational & c, const bigrational & a, const bigrational & b)
{
  if (a.is_zero())
  {
    c.num.assign(b.num);
    c.num.negate();
    c.den.assign(b.den);
  }
  else if (b.is_zero())
  {
    c.num.assign(a.num);
    c.den.assign(a.den);
  }
  else
  {
    bigint g = gcd(a.den, b.den), h;
    if (g.is_one())
    {
      // c.num   a.num * b.den - a.den * b.num;
      // ----- = -----------------------------
      // c.den         a.den * b.den;

      multiply(h, a.num, b.den);
      multiply(g, a.den, b.num);
      subtract(c.num, h, g);
      multiply(c.den, a.den, b.den);
    }
    else
    {
      // bigint t = a.num * (b.den / g) - b.num * (a.den / g);
      // bigint h = gcd(t, g);
      // c.num = t / h;
      // c.den = (a.den / g) * (b.den / h);

      bigint t, s, ss;
      divide(s, b.den, g);
      multiply(t, a.num, s);
      divide(s, a.den, g);
      multiply(ss, s, b.num);
      subtract(t, t, ss);
      h = gcd(t, g);
      divide(c.num, t, h);
      divide(t, b.den, h);
      multiply(c.den, s, t);
    }
    if (c.den.is_negative())
    {
      c.num.negate();
      c.den.negate();
    }
  }
}

void
multiply(bigrational & c, const bigrational & a, const bigrational & b)
{
  bigint g = gcd(a.num, b.den);
  bigint h = gcd(a.den, b.num);
  bigint s, t;

  switch (g.is_one() * 2 + h.is_one())
  {
    case 0:
      // c.num   (a.num / g) * (b.num / h)
      // ----- = -------------------------
      // c.den   (a.den / h) * (b.den / g)
      divide(s, a.num, g);
      divide(t, b.num, h);
      multiply(c.num, s, t);

      divide(s, a.den, h);
      divide(t, b.den, g);
      multiply(c.den, s, t);
      break;
    case 1:
      // c.num   (a.num / g) * b.num
      // ----- = -------------------
      // c.den   a.den * (b.den / g)
      divide(s, a.num, g);
      multiply(c.num, s, b.num);
      divide(t, b.den, g);
      multiply(c.den, a.den, t);
      break;
    case 2:
      // c.num   a.num * (b.num / h)
      // ----- = -------------------
      // c.den   (a.den / h) * b.den
      divide(t, b.num, h);
      multiply(c.num, a.num, t);
      divide(s, a.den, h);
      multiply(c.den, s, b.den);
      break;
    case 3:
      // c.num   a.num * b.num
      // ----- = -------------
      // c.den   a.den * b.den
      multiply(c.num, a.num, b.num);
      multiply(c.den, a.den, b.den);
      break;
  }
  if (c.den.is_negative())
  {
    c.num.negate();
    c.den.negate();
  }
}

void
divide(bigrational & c, const bigrational & a, const bigrational & b)
{
  bigint g = gcd(a.num, b.num);
  bigint h = gcd(a.den, b.den);
  bigint s, t;

  switch (g.is_one() * 2 + h.is_one())
  {
    case 0:
      // c.num   (a.num / g) * (b.den / h)
      // ----- = -------------------------
      // c.den   (a.den / h) * (b.num / g)
      divide(s, a.num, g);
      divide(t, b.den, h);
      multiply(c.num, s, t);

      divide(s, a.den, h);
      divide(t, b.num, g);
      multiply(c.den, s, t);
      break;
    case 1:
      // c.num   (a.num / g) * b.den
      // ----- = -------------------
      // c.den   a.den * (b.num / g)
      divide(s, a.num, g);
      multiply(c.num, s, b.den);
      divide(t, b.num, g);
      multiply(c.den, a.den, t);
      break;
    case 2:
      // c.num   a.num * (b.den / h)
      // ----- = -------------------
      // c.den   (a.den / h) * b.num
      divide(t, b.den, h);
      multiply(c.num, a.num, t);
      divide(s, a.den, h);
      multiply(c.den, s, b.num);
      break;
    case 3:
      // c.num   a.num * b.den
      // ----- = -------------
      // c.den   a.den * b.num
      multiply(c.num, a.num, b.den);
      multiply(c.den, a.den, b.num);
      break;
  }
  if (c.den.is_negative())
  {
    c.num.negate();
    c.den.negate();
  }
}

void
power(bigrational & c, const bigrational & a, const bigint & b)
{
  bigint n = 1, d = 1;
  if (!b.is_zero())
  {
    if (b.is_gt_zero())
    {
      power(n, a.num, b);
      power(d, a.den, b);
    }
    else
    {
      bigint abs_b = -b;
      power(n, a.den, abs_b);
      power(d, a.num, abs_b);
      if (d.is_negative())
      {
        n.negate();
        d.negate();
      }
    }
  }
  c.num.assign(n);
  c.den.assign(d);
}

void
power(bigrational & c, const bigrational & a, long i)
{
  bigint n = 1, d = 1;

  if (i)
  {
    if (i > 0)
    {
      power(n, a.num, i);
      power(d, a.den, i);
    }
    else
    {
      i = -i;
      power(n, a.den, i);
      power(d, a.num, i);
      if (d.is_negative())
      {
        n.negate();
        d.negate();
      }
    }
  }
  c.num.assign(n);
  c.den.assign(d);
}


/**
** input / output
**/

std::istream & operator >> (std::istream & in, bigrational & a)
{
  char s[10000];
  char *p = s;
  char c;

  a.num.assign_zero();
  a.den.assign_one();

  do
  {
    in.get(c);
  } while (isspace(c));
  if ((c == '+') || (c == '-'))
  {
    *p++ = c;
    do
    {
      in.get(c);
    } while (isspace(c));
  }
  if (!isdigit(c))
    std::cerr << "digit expected";
  while (isdigit(c))
  {
    *p++ = c;
    in.get(c);
  }
  *p = '\0';
  string_to_bigint(s, a.num);

  if (c != '\n')
  {
    *s = '\0';
    p = s;
    while (isspace(c))
    {
      in.get(c);
    }

    if (c == '/')
    {
      do
      {
        in.get(c);
      } while (isspace(c));
      if ((c == '+') || (c == '-'))
      {
        *p++ = c;
        do
        {
          in.get(c);
        } while (isspace(c));
      }
      if (!isdigit(c))
        std::cerr << "digit expected";
      while (isdigit(c))
      {
        *p++ = c;
        in.get(c);
      }
      *p = '\0';
      string_to_bigint(s, a.den);
    }
  }
  if (c != '\n' && c != '\r')
    in.putback(c);
  a.normalize();
  return in;
}



// Note: Remove below 4 functions since it seems we never use it in Core Library
//      Zilin Du, June 14, 2001

/**
** using fread/fwrite
**/

/*
void bigrational::
read_from_file(FILE * fp)
{
  num.read_from_file(fp);
  den.read_from_file(fp);
  this->normalize();
}

void bigrational::
write_to_file(FILE * fp)
{
  num.write_to_file(fp);
  den.write_to_file(fp);
}
*/

/**
** using fscanf/fprintf
**/

/*
void bigrational::
scan_from_file(FILE * fp)
{
  char s[10000];
  char *p = s;
  char c;

  num.assign_zero();
  den.assign_one();

  do
  {
    c = getc(fp);
  } while (isspace(c));
  if ((c == '+') || (c == '-'))
  {
    *p++ = c;
    do
    {
      c = getc(fp);
    } while (isspace(c));
  }
  else
  {
    while (isspace(c))
    {
      c = getc(fp);
    }
  }
  if (!isdigit(c))
    std::cerr << "digit expected";
  while (isdigit(c))
  {
    *p++ = c;
    c = getc(fp);
  }
  *p = '\0';
  string_to_bigint(s, num);

  if (c != '\n')
  {
    *s = '\0';
    p = s;
    while (isspace(c))
    {
      c = getc(fp);
    }

    if (c == '/')
    {
      do
      {
        c = getc(fp);
      } while (isspace(c));
      if ((c == '+') || (c == '-'))
      {
        *p++ = c;
        do
        {
          c = getc(fp);
        } while (isspace(c));
      }
      else
      {
        while (isspace(c))
        {
          c = getc(fp);
        }
      }
      if (!isdigit(c))
        std::cerr << "digit expected";
      while (isdigit(c))
      {
        *p++ = c;
        c = getc(fp);
      }
      *p = '\0';
      string_to_bigint(s, den);
    }
  }
  if (c != '\n' && c != '\r')
    ungetc(c, fp);
  this->normalize();
}

void bigrational::
print_to_file(FILE * fp)
{
  long l, k;
  char *s;

  if (num.is_zero() || den.is_one())
  {
    l = num.bit_length();
    s = new char[l / 3 + 20];
    bigint_to_string(num, s);
  }
  else
  {
    l = num.bit_length();
    k = den.bit_length();
    if (l < k)
      l = k;
    s = new char[l / 3 + 20];
    bigint_to_string(num, s);
    fputs(s, fp);
    fputs("/", fp);
    bigint_to_string(den, s);
  }
  fputs(s, fp);
  delete[] s;
}
*/


//
// non friend functions that use bigrationals
//

bool square_root(bigrational & root, const bigrational& s)
{
  if (s.is_lt_zero()) { return false; }
  bigint ns,ds;
  sqrt(ns,s.numerator());
  sqrt(ds,s.denominator());
  root = bigrational(ns,ds);
  if (root*root==s)
    return true;
  else
    return false;
}

bool cube_root(bigrational & root, const bigrational& s)
{ 
  bigint ns,ds;
  ns=s.numerator();
  if (ns.is_lt_zero()) { ns=-ns; }
  newton_root(ns,ns,3);
  newton_root(ds,s.denominator(),3);
  root = bigrational(ns,ds);
  if (s.is_lt_zero()) { root.negate(); }
  if ((root*root*root)==s)
    return true;
  else
    return false;
}

CORE_END_NAMESPACE
