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
// file          : Gmpz.h
// package       : Number_types
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Andreas Fabri
//                 Stefan Schirra
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ======================================================================
 

#ifndef CGAL_GMPZ_H
#define CGAL_GMPZ_H

#include <CGAL/basic.h>

#ifndef CGAL_CFG_NO_LOCALE
#  include <locale>
#else
#  include <cctype>
#endif // CGAL_CFG_NO_LOCALE

#include <gmp.h>

CGAL_BEGIN_NAMESPACE

class Gmpz_rep : public Rep
{
public:
  mpz_t  mpZ;

  Gmpz_rep()
  // { mpz_init_set_si(mpZ, 0); }
  { mpz_init(mpZ); }

  Gmpz_rep(mpz_t  z)
  { mpz_init_set(mpZ, z); }

  Gmpz_rep(int si)
  { mpz_init_set_si(mpZ, si); }

  Gmpz_rep(long li)
  { mpz_init_set_si(mpZ, li); }

  Gmpz_rep(unsigned long li)
  { mpz_init_set_ui(mpZ, li); }

  Gmpz_rep(double d)
  { mpz_init_set_d(mpZ, d); }

  Gmpz_rep(char* str)
  { mpz_init_set_str(mpZ, str, 10); }

  Gmpz_rep(char* str, int base)
  { mpz_init_set_str(mpZ, str, base); }

  ~Gmpz_rep()
  { mpz_clear(mpZ); }
};

class Gmpz : public Handle
{
public:
  Gmpz();

  Gmpz(const Gmpz &z);

  Gmpz(mpz_t z);

  Gmpz(int i);

  Gmpz(long l);

  Gmpz(unsigned long l);

  Gmpz(double d);

  Gmpz(char* str);
  Gmpz(char* str, int base);

  Gmpz(Gmpz_rep* R);

  ~Gmpz();

  Gmpz &operator=(const Gmpz &z);

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

  Gmpz& operator/=(const Gmpz &z);
  Gmpz operator/=(int i);

  size_t approximate_decimal_length() const;

  Gmpz_rep* ptr() const;
  double to_double() const;
  Sign sign() const;
};


inline
Gmpz_rep*
Gmpz::ptr() const
{ return static_cast<Gmpz_rep*>(PTR); }

inline   // CGAL_KERNEL_CTOR_INLINE
Gmpz::Gmpz()
{ PTR = new Gmpz_rep(0); }

inline   // CGAL_KERNEL_CTOR_INLINE
Gmpz::Gmpz(const Gmpz &z)
  : Handle(static_cast<const Handle&>(z))
{}

inline   // CGAL_KERNEL_CTOR_INLINE
Gmpz::Gmpz(mpz_t z)
{ std::cout << "OLD construction called"; PTR = new Gmpz_rep(z); }

inline   // CGAL_KERNEL_CTOR_INLINE
Gmpz::Gmpz(int i)
{ PTR = new Gmpz_rep(i); }

inline   // CGAL_KERNEL_CTOR_INLINE
Gmpz::Gmpz(long l)
{ PTR = new Gmpz_rep(l); }

inline   // CGAL_KERNEL_CTOR_INLINE
Gmpz::Gmpz(unsigned long l)
{ PTR = new Gmpz_rep(l); }

inline   // CGAL_KERNEL_CTOR_INLINE
Gmpz::Gmpz(double d)
{ PTR = new Gmpz_rep(d); }

inline   // CGAL_KERNEL_CTOR_INLINE
Gmpz::Gmpz(char* str)
{ PTR = new Gmpz_rep(str); }

inline   // CGAL_KERNEL_CTOR_INLINE
Gmpz::Gmpz(char* str, int base)
{ PTR = new Gmpz_rep(str, base); }

inline   // CGAL_KERNEL_CTOR_INLINE
Gmpz::Gmpz(Gmpz_rep* R)
{ PTR = R; }

inline
Gmpz::~Gmpz()
{}

inline
Gmpz &
Gmpz::operator=(const Gmpz &z)
{
  Handle::operator=(z);
  return *this;
}

inline
bool
Gmpz::operator==(const Gmpz &z) const
{ return ( mpz_cmp(ptr()->mpZ, z.ptr()->mpZ) == 0 ); }

inline
bool
Gmpz::operator<(const Gmpz &z) const
{ return ( mpz_cmp(ptr()->mpZ, z.ptr()->mpZ) < 0 ); }

inline
bool
Gmpz::operator<(int i) const
{ return mpz_cmp_si(ptr()->mpZ, i) < 0; }

inline
bool            /* XXX */
Gmpz::operator<=(const Gmpz &z) const
{ return ( mpz_cmp(ptr()->mpZ, z.ptr()->mpZ) <= 0 ); }


inline
bool
Gmpz::operator<=(int i) const
{ return ( mpz_cmp_si(ptr()->mpZ, i) <= 0 ); }

inline
bool             /* XXX */
Gmpz::operator>(const Gmpz &z) const
{ return ( mpz_cmp(ptr()->mpZ, z.ptr()->mpZ) > 0 ); }

inline
bool
Gmpz::operator>(int i) const
{ return ( mpz_cmp_si(ptr()->mpZ, i) > 0 ); }

inline
bool             /* XXX */
Gmpz::operator>=(const Gmpz &z) const
{ return ( mpz_cmp(ptr()->mpZ, z.ptr()->mpZ) >= 0 ); }

inline
bool
Gmpz::operator>=(int i) const
{ return ( mpz_cmp_si(ptr()->mpZ, i) >= 0 ); }

inline
bool             /* XXX */
Gmpz::operator!=(const Gmpz &z) const
{ return ! (*this == z); }

inline
bool
Gmpz::operator==(int i) const
{ return ( mpz_cmp_si(ptr()->mpZ, i) == 0 ); }

inline
bool
Gmpz::operator!=(int i) const
{ return ! (*this == i); }

inline
Gmpz
Gmpz::operator-() const
{
    Gmpz_rep* Res = new Gmpz_rep();
    mpz_neg(Res->mpZ, ptr()->mpZ);
    return Gmpz(Res);
}

inline
Gmpz
Gmpz::operator+(const Gmpz &z) const
{
    Gmpz_rep* Res = new Gmpz_rep();
    mpz_add(Res->mpZ, ptr()->mpZ, z.ptr()->mpZ);
    return Gmpz(Res);
}

inline
Gmpz
Gmpz::operator+(int i) const
{
    if (i>0)
    {
        Gmpz_rep* Res = new Gmpz_rep();
        mpz_add_ui(Res->mpZ, ptr()->mpZ, i);
        return Gmpz(Res);
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
    Gmpz_rep* Res = new Gmpz_rep();
    mpz_sub(Res->mpZ, ptr()->mpZ, z.ptr()->mpZ);
    return Gmpz(Res);
}

inline
Gmpz Gmpz::operator-(int i) const
{
    if (i>0)
    {
        Gmpz_rep* Res = new Gmpz_rep();
        mpz_sub_ui(Res->mpZ, ptr()->mpZ, i);
        return Gmpz(Res);
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
    Gmpz_rep* Res = new Gmpz_rep();
    mpz_mul(Res->mpZ, ptr()->mpZ, z.ptr()->mpZ);
    return Gmpz(Res);
}

inline
Gmpz
Gmpz::operator*(int i) const
{
    if (i>0)
    {
        Gmpz_rep* Res = new Gmpz_rep();
        mpz_mul_ui(Res->mpZ, ptr()->mpZ, i);
        return Gmpz(Res);
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
    Gmpz_rep* Res = new Gmpz_rep();
    mpz_tdiv_q(Res->mpZ, ptr()->mpZ, z.ptr()->mpZ);
    return Gmpz(Res);
}

inline
Gmpz
Gmpz::operator%(const Gmpz &z) const
{
    Gmpz_rep* Res = new Gmpz_rep();
    mpz_tdiv_r(Res->mpZ, ptr()->mpZ, z.ptr()->mpZ);
    return Gmpz(Res);
}

inline
Gmpz
Gmpz::operator/(int i) const
{
    if (i>0)
    {
        Gmpz_rep* Res = new Gmpz_rep();
        mpz_tdiv_q_ui(Res->mpZ, ptr()->mpZ, i);
        return Gmpz(Res);
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
{ return mpz_get_d(ptr()->mpZ); }

inline
io_Operator
io_tag(const Gmpz&)
{ return io_Operator(); }

inline
Sign
Gmpz::sign() const
{ return static_cast<Sign>(mpz_sgn(ptr()->mpZ)); }

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
Number_tag
number_type_tag(const Gmpz& )
{ return Number_tag(); }

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
  Gmpz_rep* Res = new Gmpz_rep();
  mpz_sqrt(Res->mpZ, z.ptr()->mpZ);
  return Gmpz(Res);
}

inline
Gmpz
gcd(const Gmpz &z1, const Gmpz &z2)
{
  Gmpz_rep* Res = new Gmpz_rep();
  mpz_gcd(Res->mpZ, z1.ptr()->mpZ, z2.ptr()->mpZ);
  return Gmpz(Res);
}

inline
Gmpz
gcd(const Gmpz &z, int i)
{
  if (i > 0)
  {
      Gmpz_rep* Res = new Gmpz_rep();
      mpz_gcd_ui(Res->mpZ, z.ptr()->mpZ, i);
      return Gmpz(Res);
  }
  return gcd(z, Gmpz(i));
}

inline
Gmpz
exact_division(const Gmpz &z1, const Gmpz &z2)
{
  Gmpz_rep* Res = new Gmpz_rep();
  mpz_divexact(Res->mpZ, z1.ptr()->mpZ, z2.ptr()->mpZ);
#ifdef CGAL_CHECK_POSTCONDITIONS
  mpz_t prod;
  mpz_init(prod);
  mpz_mul(prod, Res->mpZ, z2.ptr()->mpZ);
  CGAL_kernel_postcondition_msg(mpz_cmp(prod, z1.ptr()->mpZ) == 0,
                                "exact_division failed\n");
  mpz_clear( prod);
#endif // CGAL_CHECK_POSTCONDITIONS
  return Gmpz(Res);
}

inline
size_t
Gmpz::approximate_decimal_length() const
{ return mpz_sizeinbase(ptr()->mpZ,10); }

inline
std::ostream&
operator<<(std::ostream& os, const Gmpz &z)
{
  char *str = new char [mpz_sizeinbase(z.ptr()->mpZ,10) + 2];
  str = mpz_get_str(str, 10, z.ptr()->mpZ);
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

CGAL_END_NAMESPACE

#endif // CGAL_GMPZ_H
