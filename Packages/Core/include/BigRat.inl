/******************************************************************
 * Core Library, Version 1.5, August 2002
 * Copyright (c) 1995-2002 Exact Computation Project
 * 
 * File: BigRat.inl
 * $Id$ 
******************************************************************/
#ifdef CORE_INLINE

/**
** constructors and destructor; we could leave out some of these
**/

CORE_INLINE bigrational::bigrational()
{
  num.assign_zero();
  den.assign_one();
}

CORE_INLINE bigrational::bigrational(int i)
{
  num = i;
  den.assign_one();
}

CORE_INLINE bigrational::bigrational(long l)
{
  num = l;
  den.assign_one();
}

CORE_INLINE bigrational::bigrational(unsigned long ul)
{
  num = ul;
  den.assign_one();
}

CORE_INLINE bigrational::bigrational(const bigint & n)
{
  num.assign(n);
  den.assign_one();
}

CORE_INLINE bigrational::bigrational(const bigint & n, const bigint & d)
{
  if (d.is_zero())
  {
    lidia_error_handler("bigrational", "constructor(n,d)::division by zero.");
    assign_zero();
  }
  else
  {
   num.assign(n);
   den.assign(d);
   this->normalize();
  }
}

CORE_INLINE bigrational::bigrational(const bigrational & a)
{
  num.assign(a.num);
  den.assign(a.den);
}

CORE_INLINE bigrational::~bigrational()
{
}

/**
 ** CORE_INLINE member functions
 **/

CORE_INLINE int bigrational::
sign() const
{
  return num.sign();
}

CORE_INLINE bool bigrational::
is_positive() const
{
  if (num.sign() > 0) return true;
  return false;
}

CORE_INLINE bool bigrational::
is_negative() const
{
  if (num.sign() < 0) return true;
  return false;
}

CORE_INLINE bool bigrational::
is_zero() const
{
  if (num.sign() == 0) return true;
  return false;
}

CORE_INLINE bool bigrational::
is_gt_zero() const
{
  if (num.sign() > 0) return true;
  return false;
}

CORE_INLINE bool bigrational::
is_ge_zero() const
{
  if (num.sign() >= 0) return true;
  return false;
}

CORE_INLINE bool bigrational::
is_lt_zero() const
{
  if (num.sign() < 0) return true;
  return false;
}

CORE_INLINE bool bigrational::
is_le_zero() const
{
  if (num.sign() <= 0) return true;
  return false;
}

CORE_INLINE bool bigrational::
is_one() const
{
  if (num.is_one() && den.is_one()) return true;
  return false;
}

CORE_INLINE bool bigrational::
intify(int &i) const
{
  bigint I;
  divide(I, num, den);
  return I.intify(i);
}

CORE_INLINE bool bigrational::
longify(long &i) const
{
  bigint I;
  divide(I, num, den);
  return I.longify(i);
}

CORE_INLINE int bigrational::
abs_compare(const bigrational & a) const
{
  bigint r1, r2;
  multiply(r1, num, a.den);
  multiply(r2, den, a.num);
  return r1.abs_compare(r2);
}

CORE_INLINE int bigrational::
compare(const bigrational & a) const
{
  bigint r1, r2;
  multiply(r1, num, a.den);
  multiply(r2, den, a.num);
  return r1.compare(r2);
}

CORE_INLINE void bigrational::
normalize()
{
  int s = den.sign();
  if (s == 0)
  {
    lidia_error_handler("bigrational", "normalize()::division by zero.");
    assign_zero();
    return;
  }
  if (s < 0)
  {
    num.negate();
    den.negate();
  }
  bigint g = bgcd(num, den);
  num /= g;
  den /= g;
}

CORE_INLINE void bigrational::
absolute_value()
{
  num.absolute_value();
}

CORE_INLINE void bigrational::
negate()
{
  num.negate();
}

CORE_INLINE void bigrational::
invert()
{
  if (num.is_zero())
  {
    lidia_error_handler("bigrational", "invert()::division by zero.");
    assign_zero();
    return;
  }
  swap(num, den);
  if (den.is_negative())
  {
    num.negate();
    den.negate();
  }
}

CORE_INLINE bigint bigrational::
numerator() const
{
  return num;
}

CORE_INLINE bigint bigrational::
denominator() const
{
  return den;
}

CORE_INLINE void bigrational::
multiply_by_denominator()
{
  den.assign_one();
}

CORE_INLINE void bigrational::
assign_zero()
{
  num.assign_zero();
  den.assign_one();
}

CORE_INLINE void bigrational::
assign_one()
{
  num.assign_one();
  den.assign_one();
}

CORE_INLINE void bigrational::
assign(const bigint & n,
       const bigint & d)
{
  if (d.is_zero())
  {
    lidia_error_handler("bigrational", "assign(n,d)::division by zero.");
    assign_zero();
    return;
  }
  num.assign(n);
  den.assign(d);
  this->normalize();
}

CORE_INLINE void bigrational::
assign(const bigrational & a)
{
  num.assign(a.num);
  den.assign(a.den);
}

CORE_INLINE void bigrational::
multiply_by_2()
{
  if (den.is_even())
    den.divide_by_2();
  else
    num.multiply_by_2();
}

CORE_INLINE void bigrational::
divide_by_2()
{
  if (num.is_even())
    num.divide_by_2();
  else
    den.multiply_by_2();
}

CORE_INLINE bigint bigrational
::characteristic () const
 {
   return 0;
 }


/**
** Type checking
**/

CORE_INLINE bool
is_bigint(const bigrational & a)
{
  return a.den.is_one();
}

/**
** assignments
**/

CORE_INLINE int bigrational::operator = (int i)
{
  num = i;
  den.assign_one();
  return i;
}

CORE_INLINE long bigrational::operator = (long l)
{
  num = l;
  den.assign_one();
  return l;
}

CORE_INLINE unsigned long bigrational::operator = (unsigned long ul)
{
  num = ul;
  den.assign_one();
  return ul;
}

CORE_INLINE bigint bigrational::operator = (const bigint & a)
{
  num.assign(a);
  den.assign_one();
  return a;
}

CORE_INLINE bigrational & bigrational::operator = (const bigrational & a)
{
  if (this == &a)               // Bug fixed: avoid self-assignment
    return *this;               // Zilin Du, 07/11/01
  num.assign(a.num);
  den.assign(a.den);
  return *this;
}

/**
** comparisons
**/

CORE_INLINE bool operator == (const bigrational & a, const bigrational & b)
{ if (a.compare(b) == 0) return true; 
  return false; }

CORE_INLINE bool operator != (const bigrational & a, const bigrational & b)
{ if (a.compare(b) != 0) return true; 
  return false; }

CORE_INLINE bool operator > (const bigrational & a, const bigrational & b)
{ if (a.compare(b) > 0) return true; 
  return false; }

CORE_INLINE bool operator >= (const bigrational & a, const bigrational & b)
{ if (a.compare(b) >= 0) return true;
  return false; }

CORE_INLINE bool operator < (const bigrational & a, const bigrational & b)
{ if (a.compare(b) < 0) return true; 
  return false; }

CORE_INLINE bool operator <= (const bigrational & a, const bigrational & b)
{ if (a.compare(b) <= 0) return true;
  return false; }

/**
** operator overloading
**/

CORE_INLINE bigrational operator - (const bigrational & a)
{
  bigrational c(a);
  c.num.negate();
  return c;
}

CORE_INLINE bigrational operator + (const bigrational & a, const bigrational & b)
{
  bigrational c;
  add(c, a, b);
  return c;
}

CORE_INLINE bigrational operator - (const bigrational & a, const bigrational & b)
{
  bigrational c;
  subtract(c, a, b);
  return c;
}

CORE_INLINE bigrational operator * (const bigrational & a, const bigrational & b)
{
  bigrational c;
  multiply(c, a, b);
  return c;
}

CORE_INLINE bigrational operator / (const bigrational & a, const bigrational & b)
{
  bigrational c;
  divide(c, a, b);
  return c;
}

CORE_INLINE bigrational operator << (const bigrational & a, unsigned long ui)
{
  bigrational c;
  shift_left(c, a, ui);
  return c;
}

CORE_INLINE bigrational operator >> (const bigrational & a, unsigned long ui)
{
  bigrational c;
  shift_right(c, a, ui);
  return c;
}

CORE_INLINE bigrational & bigrational::operator += (const bigrational & a)
{
  add(*this, *this, a);
  return *this;
}

CORE_INLINE bigrational & bigrational::operator -= (const bigrational & a)
{
  subtract(*this, *this, a);
  return *this;
}

CORE_INLINE bigrational & bigrational::operator *= (const bigrational & a)
{
  multiply(*this, *this, a);
  return *this;
}

CORE_INLINE bigrational & bigrational::operator /= (const bigrational & a)
{
  divide(*this, *this, a);
  return *this;
}

CORE_INLINE bigrational & bigrational::operator <<= (unsigned long ui)
{
  shift_left(*this, *this, ui);
  return *this;
}

CORE_INLINE bigrational & bigrational::operator >>= (unsigned long ui)
{
  shift_right(*this, *this, ui);
  return *this;
}

CORE_INLINE bigrational & bigrational::operator++ ()
{
  add(num, num, den);
  return *this;
}

CORE_INLINE bigrational & bigrational::operator-- ()
{
  subtract(num, num, den);
  return *this;
}

CORE_INLINE int bigrational::operator ! ()
{
  return num.is_zero();
}

/**
 ** Procedural versions
 **/

CORE_INLINE void 
invert(bigrational & a, const bigrational & b)
{
  a.assign(b);
  a.invert();
}

CORE_INLINE void 
lidia_negate(bigrational & a, const bigrational & b)
{
  a.assign(b);
  a.negate();
}

CORE_INLINE void multiply_by_two(bigrational & c, const bigrational & a)
{
  c.assign(a);
  c.multiply_by_2();
}

CORE_INLINE void divide_by_two(bigrational & c, const bigrational & a)
{
  c.assign(a);
  c.divide_by_2();
}

CORE_INLINE void
shift_left(bigrational & c, const bigrational & a, long ui)
{
  if (ui < 0) 
  {
    lidia_error_handler("bigrational", "shift_left()::index is negative.");
    c.num.assign(a.num);
    c.den.assign(a.den);
    return;
  }
  c.den.assign(a.den);
  shift_left(c.num, a.num, ui);
  c.normalize();
}

CORE_INLINE void
shift_right(bigrational & c, const bigrational & a, long ui)
{
  if (ui < 0) 
  {
    lidia_error_handler("bigrational", "shift_right()::index is negative.");
    c.num.assign(a.num);
    c.den.assign(a.den);
    return;
  }
  c.num.assign(a.num);
  shift_left(c.den, a.den, ui);
  c.normalize();
}

CORE_INLINE void
inc(bigrational & c)
{
  add(c.num, c.num, c.den);
}

CORE_INLINE void
dec(bigrational & c)
{
  subtract(c.num, c.num, c.den);
}

CORE_INLINE void
add(bigrational & c, const bigrational & a, long i)
{
  bigint tmp;
  multiply(tmp, a.den, i);
  add(c.num, a.num, tmp);
  c.den.assign(a.den);
}

CORE_INLINE void
add(bigrational & c, long i, const bigrational & a)
{
  bigint tmp;
  multiply(tmp, a.den, i);
  add(c.num, a.num, tmp);
  c.den.assign(a.den);
}

CORE_INLINE void
subtract(bigrational & c, const bigrational & a, long i)
{
  bigint tmp;
  multiply(tmp, a.den, i);
  subtract(c.num, a.num, tmp);
  c.den.assign(a.den);
}

CORE_INLINE void
subtract(bigrational & c, long i, const bigrational & a)
{
  bigint tmp;
  multiply(tmp, a.den, i);
  subtract(c.num, a.num, tmp);
  c.num.negate();
  c.den.assign(a.den);
}

CORE_INLINE void
multiply(bigrational & c, const bigrational & a, long i)
{
  multiply(c.num, a.num, i);
  c.den.assign(a.den);
  c.normalize();
}

CORE_INLINE void
multiply(bigrational & c, long i, const bigrational & a)
{
  multiply(c.num, a.num, i);
  c.den.assign(a.den);
  c.normalize();
}

CORE_INLINE void
divide(bigrational & c, const bigrational & a, long i)
{
  c.num.assign(a.num);
  multiply(c.den, a.den, i);
  c.normalize();
}

CORE_INLINE void
divide(bigrational & c, long i, const bigrational & a)
{
  if (&c == &a)
  {
    c.invert();
    multiply(c.num, c.num, i);
  }
  else
  {
    multiply(c.num, a.den, i);
    c.den.assign(a.num);
  }
  c.normalize();
}

/**
 ** functions
 **/

CORE_INLINE bigrational
abs(const bigrational & a)
{
  bigrational c(a);
  c.absolute_value();
  return c;
}

CORE_INLINE bigrational
inverse(const bigrational & a)
{
  bigrational c(a);
  c.invert();
  return c;
}

CORE_INLINE bigint
numerator(const bigrational & a)
{
  bigint c(a.num);
  return c;
}

CORE_INLINE bigint
denominator(const bigrational & a)
{
  bigint c(a.den);
  return c;
}

CORE_INLINE bigint
round(const bigrational & a)
{
  bigrational c;
  bigint rn(a.num), rd(a.den);
  rn.multiply_by_2();
  add(rn, rn, rd);
  rd.multiply_by_2();
  c.num.assign(rn);
  c.den.assign(rd);
  return floor(c);
}

CORE_INLINE bigint
floor(const bigrational & a)
{
  bigint q, r;
  div_rem(q, r, abs(a.num), abs(a.den));
  if (a.sign() < 0 && r.sign() != 0)
  {
    if (a.sign() < 0)
      q.negate();
    dec(q);
  }
  else
  {
    if (a.sign() < 0)
      q.negate();
  }
  return q;
}

CORE_INLINE bigint
ceiling(const bigrational & a)
{
  bigint q, r;
  div_rem(q, r, abs(a.num), abs(a.den));
  if (a.sign() >= 0 && r.sign() != 0)
    inc(q);
  if (a.sign() < 0)
    q.negate();
  return q;
}

CORE_INLINE bigint
truncate(const bigrational & a)
{
  bigint q, r;
  div_rem(q, r, abs(a.num), abs(a.den));
  if (a.sign() < 0)
    q.negate();
  return q;
}

CORE_INLINE double
dbl(const bigrational & a)
{
  if (a.num.is_zero())
    return 0.0;

  long ln = a.num.bit_length();
  long ld = a.den.bit_length();
  long en = 0, ed = 0;

  if (ln > 1023)
    en = ln - 1023;
  if (ld > 1023)
    ed = ld - 1023;

  bigint an = a.num >> en;
  bigint ad = a.den >> ed;

  // an = (a/2^en)*2^en
  // ad = (a/2^ed)*2^ed

  double d = dbl(an) / dbl(ad);
  d = ldexp(d, (int)(en - ed));

  return d;
}

CORE_INLINE void
square(bigrational & a, const bigrational & b)
{
  square(a.num, b.num);
  square(a.den, b.den);
}

CORE_INLINE void
swap(bigrational & a, bigrational & b)
{
  swap(a.num, b.num);
  swap(a.den, b.den);
}

CORE_INLINE std::ostream & operator << (std::ostream & out, const bigrational & a)
{
  bigint n = a.num;
  bigint d = a.den;

  if (n.is_zero() || d.is_one())
    out << n;
  else
    out << n << "/" << d;
  return out;
}

CORE_INLINE int
string_to_bigrational(char *s, char *t, bigrational & a)
{
  long i = string_to_bigint(s, a.num);
  long j = string_to_bigint(t, a.den);
  a.normalize();
  return (i + j);
}

CORE_INLINE int
bigrational_to_string(const bigrational & a, char *s, char *t)
{
  int i = bigint_to_string(a.num, s);
  int j = bigint_to_string(a.den, t);
  return (i + j);
}

#endif

