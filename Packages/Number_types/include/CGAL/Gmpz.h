// Copyright (c) 1999,2003  Utrecht University (The Netherlands),
// ETH Zurich (Switzerland), Freie Universitaet Berlin (Germany),
// INRIA Sophia-Antipolis (France), Martin-Luther-University Halle-Wittenberg
// (Germany), Max-Planck-Institute Saarbruecken (Germany), RISC Linz (Austria),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Andreas Fabri, Stefan Schirra, Sylvain Pion
 

#ifndef CGAL_GMPZ_H
#define CGAL_GMPZ_H

#include <CGAL/basic.h>
#include <CGAL/Handle_for.h>
#include <CGAL/Quotient.h>
#include <CGAL/double.h> 
#include <CGAL/Interval_arithmetic.h>

#include <string>
#ifndef CGAL_CFG_NO_LOCALE
#  include <locale>
#else
#  include <cctype>
#endif

#include <gmp.h>

CGAL_BEGIN_NAMESPACE

class Gmpz_rep
{
public:
  mpz_t  mpZ;

  Gmpz_rep()
  { mpz_init(mpZ); }

  Gmpz_rep(const mpz_t z)
  { mpz_init_set(mpZ, z); }

  Gmpz_rep(const Gmpz_rep & g)
  { mpz_init_set(mpZ, g.mpZ); }

  Gmpz_rep & operator= (const Gmpz_rep & g)
  {
      if (&g != this) {
	  mpz_clear(mpZ);
	  mpz_init_set(mpZ, g.mpZ);
      }
      return *this;
  }

  Gmpz_rep(int si)
  { mpz_init_set_si(mpZ, si); }

  Gmpz_rep(long li)
  { mpz_init_set_si(mpZ, li); }

  Gmpz_rep(unsigned long li)
  { mpz_init_set_ui(mpZ, li); }

  Gmpz_rep(double d)
  { mpz_init_set_d(mpZ, d); }

  Gmpz_rep(const char * const str)
  { mpz_init_set_str(mpZ, str, 10); }

  Gmpz_rep(const char * const str, int base)
  { mpz_init_set_str(mpZ, str, base); }

  ~Gmpz_rep()
  { mpz_clear(mpZ); }
};

class Gmpz
  : public Handle_for<Gmpz_rep>
{
  typedef Handle_for<Gmpz_rep> Base;
public:
  typedef Tag_true  Has_gcd;
  typedef Tag_false Has_division;
  typedef Tag_true  Has_sqrt;

  Gmpz() // {} we can't do that since the non-const mpz() is called.
    : Base(Gmpz_rep()) {}

  Gmpz(const mpz_t z)
    : Base(Gmpz_rep(z)) {}

  Gmpz(int i)
    : Base(Gmpz_rep(i)) {}

  Gmpz(long l)
    : Base(Gmpz_rep(l)) {}

  Gmpz(unsigned long l)
    : Base(Gmpz_rep(l)) {}

  Gmpz(double d)
    : Base(Gmpz_rep(d)) {}

  Gmpz(const std::string& str)
    : Base(Gmpz_rep(str.c_str())) {}

  Gmpz(const std::string& str, int base)
    : Base(Gmpz_rep(str.c_str(), base)) {}

  Gmpz operator-() const;

  Gmpz operator%(const Gmpz &z) const;

  Gmpz& operator+=(const Gmpz &z);
  Gmpz operator+=(int i);

  Gmpz& operator-=(const Gmpz &z);
  Gmpz operator-=(int i);

  Gmpz& operator*=(const Gmpz &z);
  Gmpz operator*=(int i);

  Gmpz& operator%=(const Gmpz &z);

  Gmpz& operator/=(const Gmpz &z);
  Gmpz operator/=(int i);

  size_t approximate_decimal_length() const;

  double to_double() const;
  Sign sign() const;

  const mpz_t & mpz() const { return Ptr()->mpZ; }
  mpz_t & mpz() { return ptr()->mpZ; }
};


inline
bool
operator==(const Gmpz &a, const Gmpz &b)
{ return mpz_cmp(a.mpz(), b.mpz()) == 0; }

inline
bool
operator<(const Gmpz &a, const Gmpz &b)
{ return mpz_cmp(a.mpz(), b.mpz()) < 0; }

inline
bool
operator<=(const Gmpz &a, const Gmpz &b)
{ return ! (b < a); }

inline
bool
operator>(const Gmpz &a, const Gmpz &b)
{ return b < a; }

inline
bool
operator>=(const Gmpz &a, const Gmpz &b)
{ return ! (a < b); }

inline
bool
operator!=(const Gmpz &a, const Gmpz &b)
{ return ! (a == b); }

// mixed operators.
inline
bool
operator<(const Gmpz &a, int b)
{ return mpz_cmp_si(a.mpz(), b) < 0; }

inline
bool
operator<(int a, const Gmpz &b)
{ return mpz_cmp_si(b.mpz(), a) > 0; }

inline
bool
operator==(const Gmpz &a, int b)
{ return mpz_cmp_si(a.mpz(), b) == 0; }

inline
bool
operator==(int a, const Gmpz &b)
{ return b == a; }

inline
bool
operator<=(const Gmpz &a, int b)
{ return ! (b < a); }

inline
bool
operator<=(int a, const Gmpz &b)
{ return ! (b < a); }

inline
bool
operator>(const Gmpz &a, int b)
{ return b < a; }

inline
bool
operator>(int a, const Gmpz &b)
{ return b < a; }

inline
bool
operator>=(const Gmpz &a, int b)
{ return ! (a < b); }

inline
bool
operator>=(int a, const Gmpz &b)
{ return ! (a < b); }

inline
bool
operator!=(const Gmpz &a, int b)
{ return ! (a == b); }

inline
bool
operator!=(int a, const Gmpz &b)
{ return ! (a == b); }


inline
Gmpz
Gmpz::operator-() const
{
    Gmpz Res;
    mpz_neg(Res.mpz(), mpz());
    return Res;
}


inline
Gmpz
operator+(const Gmpz &a, const Gmpz &b)
{
    Gmpz Res;
    mpz_add(Res.mpz(), a.mpz(), b.mpz());
    return Res;
}

inline
Gmpz
operator+(const Gmpz &a, int b)
{
    if (b>0)
    {
        Gmpz Res;
        mpz_add_ui(Res.mpz(), a.mpz(), b);
        return Res;
    }
    return a + Gmpz(b);
}

inline
Gmpz&
Gmpz::operator+=(const Gmpz &z)
{
    *this = *this + z;
    return *this;
}

inline
Gmpz
Gmpz::operator+=(int i)
{
    *this = *this + Gmpz(i);
    return *this;
}

inline
Gmpz
operator-(const Gmpz &a, const Gmpz &b)
{
    Gmpz Res;
    mpz_sub(Res.mpz(), a.mpz(), b.mpz());
    return Res;
}

inline
Gmpz
operator-(const Gmpz &a, int b)
{
    if (b>0)
    {
        Gmpz Res;
        mpz_sub_ui(Res.mpz(), a.mpz(), b);
        return Res;
    }
    return a - Gmpz(b);
}

inline
Gmpz&
Gmpz::operator-=(const Gmpz &z)
{
    *this = *this - z;
    return *this;
}

inline
Gmpz
Gmpz::operator-=(int i)
{
    *this = *this - Gmpz(i);
    return *this;
}

inline
Gmpz
operator*(const Gmpz &a, const Gmpz &b)
{
    Gmpz Res;
    mpz_mul(Res.mpz(), a.mpz(), b.mpz());
    return Res;
}

inline
Gmpz
operator*(const Gmpz &a, int b)
{
    if (b>0)
    {
        Gmpz Res;
        mpz_mul_ui(Res.mpz(), a.mpz(), b);
        return Res;
    }
    return a * Gmpz(b);
}

inline
Gmpz&
Gmpz::operator*=(const Gmpz &z)
{
    *this = *this * z;
    return *this;
}

inline
Gmpz
Gmpz::operator*=(int i)
{
    *this = *this * Gmpz(i);
    return *this;
}

inline
Gmpz
operator/(const Gmpz &a, const Gmpz &b)
{
    CGAL_precondition(b != 0);
    Gmpz Res;
    mpz_tdiv_q(Res.mpz(), a.mpz(), b.mpz());
    return Res;
}

inline
Gmpz
operator/(const Gmpz &a, int b)
{
    if (b>0)
    {
        Gmpz Res;
        mpz_tdiv_q_ui(Res.mpz(), a.mpz(), b);
        return Res;
    }
    return a / Gmpz(b);
}

inline
Gmpz&
Gmpz::operator/=(const Gmpz &z)
{
    *this = *this / z;
    return *this;
}

inline
Gmpz
Gmpz::operator/=(int i)
{
    *this = *this / Gmpz(i);
    return *this;
}

inline
Gmpz
operator+(int i, const Gmpz &z)
{ return z + i; }

inline
Gmpz
operator-(int i, const Gmpz &z)
{ return Gmpz(i) - z; }

inline
Gmpz
operator*(int i, const Gmpz &z)
{ return z * i; }

inline
Gmpz
operator/(int i, const Gmpz &z)
{ return Gmpz(i) / z; }

inline
Gmpz
Gmpz::operator%(const Gmpz &z) const
{
    Gmpz Res;
    mpz_tdiv_r(Res.mpz(), mpz(), z.mpz());
    return Res;
}

inline
Gmpz&
Gmpz::operator%=(const Gmpz &z)
{
    *this = *this % z;
    return *this;
}

inline
double
Gmpz::to_double() const
{ return mpz_get_d(mpz()); }

inline
io_Operator
io_tag(const Gmpz&)
{ return io_Operator(); }

inline
Sign
Gmpz::sign() const
{ return static_cast<Sign>(mpz_sgn(mpz())); }

inline
double
to_double(const Gmpz &z)
{ return z.to_double(); }

inline
Sign
sign(const Gmpz &z)
{ return z.sign(); }

inline
bool
is_valid(const Gmpz &)
{ return true; }

inline
bool
is_finite(const Gmpz &)
{ return true; }

inline
Gmpz
sqrt(const Gmpz &z)
{
  Gmpz Res;
  mpz_sqrt(Res.mpz(), z.mpz());
  return Res;
}

inline
Gmpz
div(const Gmpz &z1, const Gmpz &z2)
{
  return z1 / z2;
}

inline
Gmpz
gcd(const Gmpz &z1, const Gmpz &z2)
{
  Gmpz Res;
  mpz_gcd(Res.mpz(), z1.mpz(), z2.mpz());
  return Res;
}

inline
Gmpz
gcd(const Gmpz &z, int i)
{
  if (i > 0)
  {
      Gmpz Res;
      mpz_gcd_ui(Res.mpz(), z.mpz(), i);
      return Res;
  }
  return gcd(z, Gmpz(i));
}

inline
Gmpz
exact_division(const Gmpz &z1, const Gmpz &z2)
{
  Gmpz Res;
  mpz_divexact(Res.mpz(), z1.mpz(), z2.mpz());
#ifdef CGAL_CHECK_POSTCONDITIONS
  mpz_t prod;
  mpz_init(prod);
  mpz_mul(prod, Res.mpz(), z2.mpz());
  CGAL_kernel_postcondition_msg(mpz_cmp(prod, z1.mpz()) == 0,
                                "exact_division failed\n");
  mpz_clear( prod);
#endif // CGAL_CHECK_POSTCONDITIONS
  return Res;
}

inline
size_t
Gmpz::approximate_decimal_length() const
{ return mpz_sizeinbase(mpz(),10); }

inline
std::ostream&
operator<<(std::ostream& os, const Gmpz &z)
{
  char *str = new char [mpz_sizeinbase(z.mpz(),10) + 2];
  str = mpz_get_str(str, 10, z.mpz());
  os << str ;
  delete[] str;
  return os;
}

inline
std::istream&
operator>>(std::istream& is, Gmpz &z)
{
  int negative = 0;
  const int null = '0';
  char c;

#ifndef CGAL_CFG_NO_LOCALE
  while (is.get(c) && std::isspace(c, std::locale::classic() ))
#else
  while (is.get(c) && CGAL_CLIB_STD::isspace(c))
#endif // CGAL_CFG_NO_LOCALE
  {}

  if (c == '-')
  {
        negative = 1;
#ifndef CGAL_CFG_NO_LOCALE
        while (is.get(c) && std::isspace(c, std::locale::classic() ))
#else
        while (is.get(c) && CGAL_CLIB_STD::isspace(c))
#endif // CGAL_CFG_NO_LOCALE
        {}
  }
#ifndef CGAL_CFG_NO_LOCALE
  if (std::isdigit(c, std::locale::classic() ))
#else
  if (std::isdigit(c))
#endif // CGAL_CFG_NO_LOCALE
  {
        z = c - null;
#ifndef CGAL_CFG_NO_LOCALE
        while (is.get(c) && std::isdigit(c, std::locale::classic() ))
#else
        while (is.get(c) && std::isdigit(c))
#endif // CGAL_CFG_NO_LOCALE
        {
            z = 10*z + (c-null);
        }
  }
  if (is)
  {
        is.putback(c);
  }
  if (sign(z) != static_cast<Sign>(0) && negative)
  {
        z = -z;
  }
  return is;
}

inline
std::pair<double,double>
to_interval (const Gmpz & z)
{
  // GMP returns the closest double (seen in the code).
  Protect_FPU_rounding<true> P(CGAL_FE_TONEAREST);
  double app = CGAL::to_double(z);
  // If it's lower than 2^53, then it's exact.
  if (CGAL_CLIB_STD::fabs(app) < double(1<<26)*double(1<<27))
      return to_interval(app);
  FPU_set_cw(CGAL_FE_UPWARD);
  Interval_nt<false> approx(app);
  approx += Interval_nt<false>::smallest();
  return approx.pair();
}

inline
double to_double(const Quotient<Gmpz>& quot)
{
  mpq_t  mpQ;
  mpq_init(mpQ);
  const Gmpz& n = quot.numerator();
  const Gmpz& d = quot.denominator();
  mpz_set(mpq_numref(mpQ), n.mpz());
  mpz_set(mpq_denref(mpQ), d.mpz());
    
  mpq_canonicalize(mpQ);
  
  double ret = mpq_get_d(mpQ);
  mpq_clear(mpQ);
  return ret;
}

CGAL_END_NAMESPACE

#endif // CGAL_GMPZ_H
