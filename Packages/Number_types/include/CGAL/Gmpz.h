// ======================================================================
//
// Copyright (c) 1999 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
// 
// release       : 
// release_date  : 
// 
// file          : include/CGAL/Gmpz.h
// package       : Number_types
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Andreas Fabri, Stefan Schirra, Sylvain Pion
//
// coordinator   : MPI, Saarbruecken
// ======================================================================
 

#ifndef CGAL_GMPZ_H
#define CGAL_GMPZ_H

#include <CGAL/basic.h>
#include <CGAL/Handle_for.h>

#ifndef CGAL_CFG_NO_LOCALE
#  include <locale>
#else
#  include <cctype>
#endif // CGAL_CFG_NO_LOCALE

#include <gmp.h>

CGAL_BEGIN_NAMESPACE

class Gmpz_rep
{
public:
  mpz_t  mpZ;

  Gmpz_rep()
  { mpz_init(mpZ); }

  Gmpz_rep(mpz_t z)
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

  Gmpz()
    : Base(Gmpz_rep()) {}

  Gmpz(mpz_t z)
    : Base(Gmpz_rep(z)) {}

  Gmpz(int i)
    : Base(Gmpz_rep(i)) {}

  Gmpz(long l)
    : Base(Gmpz_rep(l)) {}

  Gmpz(unsigned long l)
    : Base(Gmpz_rep(l)) {}

  Gmpz(double d)
    : Base(Gmpz_rep(d)) {}

  Gmpz(const char* const str)
    : Base(Gmpz_rep(str)) {}

  Gmpz(const char* const str, int base)
    : Base(Gmpz_rep(str, base)) {}

  bool operator==(const Gmpz &z) const;
  bool operator==(int i) const;

  bool operator!=(const Gmpz &z) const;
  bool operator!=(int i) const;

  bool operator<(const Gmpz &z) const;
  bool operator<(int i) const;

  bool operator<=(const Gmpz &z) const;
  bool operator<=(int i) const;

  bool operator>(const Gmpz &z) const;
  bool operator>(int i) const;

  bool operator>=(const Gmpz &z) const;
  bool operator>=(int i) const;

  Gmpz operator-() const;

  Gmpz operator+(const Gmpz &z) const;
  Gmpz operator+(int i) const;

  Gmpz operator-(const Gmpz &z) const;
  Gmpz operator-(int i) const;

  Gmpz operator*(const Gmpz &z) const;
  Gmpz operator*(int i) const;

  Gmpz operator%(const Gmpz &z) const;

  Gmpz operator/(const Gmpz &z) const;
  Gmpz operator/(int i) const;

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
Gmpz::operator==(const Gmpz &z) const
{ return mpz_cmp(mpz(), z.mpz()) == 0; }

inline
bool
Gmpz::operator<(const Gmpz &z) const
{ return mpz_cmp(mpz(), z.mpz()) < 0; }

inline
bool
Gmpz::operator<(int i) const
{ return mpz_cmp_si(mpz(), i) < 0; }

inline
bool
Gmpz::operator<=(const Gmpz &z) const
{ return mpz_cmp(mpz(), z.mpz()) <= 0; }

inline
bool
Gmpz::operator<=(int i) const
{ return mpz_cmp_si(mpz(), i) <= 0; }

inline
bool
Gmpz::operator>(const Gmpz &z) const
{ return mpz_cmp(mpz(), z.mpz()) > 0; }

inline
bool
Gmpz::operator>(int i) const
{ return mpz_cmp_si(mpz(), i) > 0; }

inline
bool
Gmpz::operator>=(const Gmpz &z) const
{ return mpz_cmp(mpz(), z.mpz()) >= 0; }

inline
bool
Gmpz::operator>=(int i) const
{ return mpz_cmp_si(mpz(), i) >= 0; }

inline
bool
Gmpz::operator!=(const Gmpz &z) const
{ return ! (*this == z); }

inline
bool
Gmpz::operator==(int i) const
{ return mpz_cmp_si(mpz(), i) == 0; }

inline
bool
Gmpz::operator!=(int i) const
{ return ! (*this == i); }

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
Gmpz::operator+(const Gmpz &z) const
{
    Gmpz Res;
    mpz_add(Res.mpz(), mpz(), z.mpz());
    return Res;
}

inline
Gmpz
Gmpz::operator+(int i) const
{
    if (i>0)
    {
        Gmpz Res;
        mpz_add_ui(Res.mpz(), mpz(), i);
        return Res;
    }
    return *this + Gmpz(i);
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
Gmpz::operator-(const Gmpz &z) const
{
    Gmpz Res;
    mpz_sub(Res.mpz(), mpz(), z.mpz());
    return Res;
}

inline
Gmpz Gmpz::operator-(int i) const
{
    if (i>0)
    {
        Gmpz Res;
        mpz_sub_ui(Res.mpz(), mpz(), i);
        return Res;
    }
    return *this - Gmpz(i);
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
Gmpz::operator*(const Gmpz &z) const
{
    Gmpz Res;
    mpz_mul(Res.mpz(), mpz(), z.mpz());
    return Res;
}

inline
Gmpz
Gmpz::operator*(int i) const
{
    if (i>0)
    {
        Gmpz Res;
        mpz_mul_ui(Res.mpz(), mpz(), i);
        return Res;
    }
    return *this * Gmpz(i);
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
Gmpz::operator/(const Gmpz &z) const
{
    Gmpz Res;
    mpz_tdiv_q(Res.mpz(), mpz(), z.mpz());
    return Res;
}

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
Gmpz
Gmpz::operator/(int i) const
{
    if (i>0)
    {
        Gmpz Res;
        mpz_tdiv_q_ui(Res.mpz(), mpz(), i);
        return Res;
    }
    return *this / Gmpz(i);
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
        z = c - '0';
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
Interval_base
to_interval (const Gmpz & z)
{
  // GMP returns the closest double (seen in the code).
  Protect_FPU_rounding<true> P(CGAL_FE_TONEAREST);
  double app = CGAL::to_double(z);
  // If it's lower than 2^53, then it's exact.
  if (CGAL_CLIB_STD::fabs(app) < double(1<<26)*double(1<<27))
      return app;
  FPU_set_cw(CGAL_FE_UPWARD);
  return Interval_nt_advanced(app) + Interval_base::Smallest;
}

namespace NTS {
  inline
  Gmpz
  gcd( const Gmpz& n1, const Gmpz& n2)
  { 
    return CGAL::gcd(n1, n2);
  }
}

CGAL_END_NAMESPACE

#endif // CGAL_GMPZ_H
