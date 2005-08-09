// Copyright (c) 2002,2003  Utrecht University (The Netherlands),
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
// Author(s)     : Andreas Fabri, Sylvain Pion


#ifndef CGAL_GMPQ_H
#define CGAL_GMPQ_H

#include <CGAL/basic.h>
#include <CGAL/Handle_for.h>
#include <CGAL/Gmpz.h>
#include <CGAL/Interval_nt.h>

#include <utility>
#include <string>

#include <gmp.h>
#include <mpfr.h>

#include <boost/operators.hpp>

CGAL_BEGIN_NAMESPACE

// TODO : add mixed operators with Gmpz.

// Wrapper around mpq_t to get the destructor call mpq_clear.
// Contrary to mpz_t, there are no mpq_init_set_* functions,
// so we simply call mpq_init() here.
struct Gmpq_rep
{
  mpq_t mpQ;

  Gmpq_rep()  { mpq_init(mpQ); }
  ~Gmpq_rep() { mpq_clear(mpQ); }

private:
  // Make sure it does not get accidentally copied.
  Gmpq_rep(const Gmpq_rep &);
  Gmpq_rep & operator= (const Gmpq_rep &);
};


class Gmpq
  : Handle_for<Gmpq_rep>,
    boost::ordered_field_operators1< Gmpq
  , boost::ordered_field_operators2< Gmpq, int
    > >
{
  typedef Handle_for<Gmpq_rep> Base;
public:
  typedef Tag_false  Has_gcd;
  typedef Tag_true   Has_division;
  typedef Tag_false  Has_sqrt;

  typedef Tag_true   Has_exact_ring_operations;
  typedef Tag_true   Has_exact_division;
  typedef Tag_false  Has_exact_sqrt;

  Gmpq() {}

  Gmpq(const mpq_t q)
  { mpq_set(mpq(), q); }

  Gmpq(int n)
  { mpq_set_si(mpq(), n, 1); }

  Gmpq(long n)
  { mpq_set_si(mpq(), n, 1); }

  Gmpq(unsigned long n)
  { mpq_set_ui(mpq(), n, 1); }

  Gmpq(const Gmpz& n)
  { mpq_set_z(mpq(), n.mpz()); }

  Gmpq(int n, int d)
  {
    if (d < 0) {
      n = -n;
      d = -d;
    }
    mpq_set_si(mpq(), n, d);
    mpq_canonicalize(mpq());
  }

  Gmpq(signed long n, unsigned long d)
  {
    mpq_set_si(mpq(), n, d);
    mpq_canonicalize(mpq());
  }

  Gmpq(unsigned long n, unsigned long d)
  {
    mpq_set_ui(mpq(), n, d);
    mpq_canonicalize(mpq());
  }

  Gmpq(const Gmpz& n, const Gmpz& d)
  {
    mpz_set(mpq_numref(mpq()), n.mpz());
    mpz_set(mpq_denref(mpq()), d.mpz());
    mpq_canonicalize(mpq());
  }

  Gmpq(double d)
  { mpq_set_d(mpq(), d); }

  Gmpq(const std::string& str, int base = 10)
  {
    mpq_set_str(mpq(), str.c_str(), base);
    mpq_canonicalize(mpq());
  }


  Gmpz numerator() const
  { return Gmpz(mpq_numref(mpq())); }

  Gmpz denominator() const
  { return Gmpz(mpq_denref(mpq())); }

  Gmpq operator-() const;

  Gmpq& operator+=(const Gmpq &z);
  Gmpq& operator-=(const Gmpq &z);
  Gmpq& operator*=(const Gmpq &z);
  Gmpq& operator/=(const Gmpq &z);

  double to_double() const;
  Sign sign() const;

  const mpq_t & mpq() const { return Ptr()->mpQ; }
  mpq_t & mpq() { return ptr()->mpQ; }
};


inline
bool
operator==(const Gmpq &a, const Gmpq &b)
{ return mpq_equal(a.mpq(), b.mpq()); }

inline
bool
operator<(const Gmpq &a, const Gmpq &b)
{ return mpq_cmp(a.mpq(), b.mpq()) < 0; }


// mixed operators.
inline
bool
operator<(const Gmpq &a, int b)
{ return mpq_cmp_si(a.mpq(), b, 1) < 0; }

inline
bool
operator>(const Gmpq &a, int b)
{ return mpq_cmp_si(a.mpq(), b, 1) > 0; }

inline
bool
operator==(const Gmpq &a, int b)
{ return mpq_cmp_si(a.mpq(), b, 1) == 0; }


inline
Gmpq
Gmpq::operator-() const
{
    Gmpq Res;
    mpq_neg(Res.mpq(), mpq());
    return Res;
}


inline
Gmpq&
Gmpq::operator+=(const Gmpq &z)
{
    Gmpq Res;
    mpq_add(Res.mpq(), mpq(), z.mpq());
    swap(Res);
    return *this;
}

inline
Gmpq&
Gmpq::operator-=(const Gmpq &z)
{
    Gmpq Res;
    mpq_sub(Res.mpq(), mpq(), z.mpq());
    swap(Res);
    return *this;
}

inline
Gmpq&
Gmpq::operator*=(const Gmpq &z)
{
    Gmpq Res;
    mpq_mul(Res.mpq(), mpq(), z.mpq());
    swap(Res);
    return *this;
}

inline
Gmpq&
Gmpq::operator/=(const Gmpq &z)
{
    CGAL_precondition(z != 0);
    Gmpq Res;
    mpq_div(Res.mpq(), mpq(), z.mpq());
    swap(Res);
    return *this;
}


inline
double
Gmpq::to_double() const
{ return mpq_get_d(mpq()); }

inline
io_Operator
io_tag(const Gmpq&)
{ return io_Operator(); }

inline
Sign
Gmpq::sign() const
{ return static_cast<Sign>(mpq_sgn(mpq())); }


inline
double
to_double(const Gmpq &z)
{ return z.to_double(); }

inline
Sign
sign(const Gmpq &z)
{ return z.sign(); }

inline
bool
is_valid(const Gmpq &)
{ return true; }

inline
bool
is_finite(const Gmpq &)
{ return true; }



inline
std::ostream&
operator<<(std::ostream& os, const Gmpq &z)
{
  os << z.numerator() << "/" << z.denominator();
  return os;
}

inline
std::istream&
operator>>(std::istream& is, Gmpq &z)
{
  char c;
  Gmpz n, d;
  is >> n;
  is >> c;
  CGAL_assertion(!is || c == '/');
  is >> d;
  if (is)
    z = Gmpq(n,d);

  return is;
}

inline
std::pair<double, double>
to_interval (const Gmpq& z)
{
    mpfr_t x;
    mpfr_init2 (x, 53); /* Assume IEEE-754 */
    mpfr_set_q (x, z.mpq(), GMP_RNDD);
    double i = mpfr_get_d (x, GMP_RNDD); /* EXACT but can overflow */
    mpfr_set_q (x, z.mpq(), GMP_RNDU);
    double s = mpfr_get_d (x, GMP_RNDU); /* EXACT but can overflow */
    mpfr_clear (x);
    return std::pair<double, double>(i, s);
}

template <>
struct Rational_traits<Gmpq> {
  typedef Gmpz RT;
  RT   numerator     (const Gmpq & r) const { return r.numerator(); }
  RT   denominator   (const Gmpq & r) const { return r.denominator(); }
  Gmpq make_rational (const RT & n, const RT & d) const
  { return Gmpq(n, d); }
  Gmpq make_rational (const Gmpq & n, const Gmpq & d) const
  { return n / d; }
};

CGAL_END_NAMESPACE

#endif // CGAL_GMPQ_H
