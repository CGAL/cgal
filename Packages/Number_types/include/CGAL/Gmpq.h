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
#include <CGAL/Interval_arithmetic.h>

#include <utility>
#include <string>
#include <gmp.h>
#include <mpfr.h>


CGAL_BEGIN_NAMESPACE

class Gmpq_rep
{
public:

  mpq_t  mpQ;

  Gmpq_rep()
  { 
    mpq_init(mpQ); 
  }

  Gmpq_rep(const mpq_t z)
  { 
    mpq_init(mpQ); 
    mpq_set(mpQ, z);
  }

  Gmpq_rep(const Gmpq_rep & g)
  { 
    mpq_init(mpQ); 
    mpq_set(mpQ, g.mpQ);
  }

  Gmpq_rep & operator= (const Gmpq_rep & g)
  {
      if (&g != this) {
	  mpq_clear(mpQ);
	  mpq_set(mpQ, g.mpQ);
      }
      return *this;
  }

  Gmpq_rep(int si)
  { 
    mpq_init(mpQ); 
    mpq_set_si(mpQ, si, 1); 
  }

  Gmpq_rep(long si)
  { 
    mpq_init(mpQ); 
    mpq_set_si(mpQ, si, 1); 
  }


  Gmpq_rep(unsigned long ui)
  { 
    mpq_init(mpQ); 
    mpq_set_ui(mpQ, ui, 1); 
  }

  Gmpq_rep(const Gmpz& z)
  { 
    mpq_init(mpQ); 
    mpq_set_z(mpQ, z.mpz()); 
  }

  Gmpq_rep(unsigned long int ui1, unsigned long int ui2)
  { 
    mpq_init(mpQ); 
    mpq_set_ui(mpQ, ui1, ui2); 
    mpq_canonicalize(mpQ);
  }
  
  Gmpq_rep(signed long int si, unsigned long int ui)
  { 
    mpq_init(mpQ); 
    mpq_set_si(mpQ, si, ui);
    mpq_canonicalize(mpQ);
  }
  
  Gmpq_rep(int num, int den)
  { 
    mpq_init(mpQ); 
    if(den < 0) {
      num = -num;
      den = -den;
    }
    mpq_set_si(mpQ, num, den);
    mpq_canonicalize(mpQ);
  }

  Gmpq_rep(const Gmpz& n, const Gmpz& d)
  { 
    mpq_init(mpQ); 
    mpz_set(mpq_numref(mpQ), n.mpz());
    mpz_set(mpq_denref(mpQ), d.mpz());
    
    mpq_canonicalize(mpQ);
  }

  Gmpq_rep(double d)
  { 
    mpq_init(mpQ);
    mpq_set_d(mpQ, d); 
  }

  Gmpq_rep(const std::string& str)
  { 
    mpq_init(mpQ);
    mpq_set_str(mpQ, str.c_str(), 10);
    mpq_canonicalize(mpQ); 
  }

  Gmpq_rep(const std::string& str, int base)
  { 
    mpq_init(mpQ);
    mpq_set_str(mpQ, str.c_str(), base);
    mpq_canonicalize(mpQ);
  }

  ~Gmpq_rep()
  { mpq_clear(mpQ); }
};


class Gmpq
  : public Handle_for<Gmpq_rep>
{
  typedef Handle_for<Gmpq_rep> Base;
public:
  typedef Tag_false  Has_gcd;
  typedef Tag_true Has_division;
  typedef Tag_false  Has_sqrt;


  Gmpq() // {} we can't do that since the non-const mpq() is called.
    : Base(Gmpq_rep()) {}

  Gmpq(const mpq_t z)
    : Base(Gmpq_rep(z)) {}

  Gmpq(int n)
    : Base(Gmpq_rep(n)) {}

  Gmpq(long n)
    : Base(Gmpq_rep(n)) {}

  Gmpq(unsigned long n)
    : Base(Gmpq_rep(n)) {}

  Gmpq(const Gmpz& n)
    : Base(Gmpq_rep(n)) {}

  Gmpq(int n, int d)
    : Base(Gmpq_rep(n, d)) {}
  
  Gmpq(signed long n, unsigned long d)
    : Base(Gmpq_rep(n, d)) {}

  Gmpq(unsigned long n, unsigned long d)
    : Base(Gmpq_rep(n, d)) {}

  Gmpq(const Gmpz& n, const Gmpz& d)
    : Base(Gmpq_rep(n,d)) {}

  Gmpq(double d)
    : Base(Gmpq_rep(d)) {}
  
  Gmpq(const std::string& str)
    : Base(Gmpq_rep(str)) {}
  
   Gmpq(const std::string& str, int base)
    : Base(Gmpq_rep(str, base)) {}


  Gmpz numerator() const
  {
    return Gmpz(mpq_numref(mpq()));
  }

  Gmpz denominator() const
  {
    return Gmpz(mpq_denref(mpq()));
    
  }

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

inline
bool
operator<=(const Gmpq &a, const Gmpq &b)
{ return ! (b < a); }

inline
bool
operator>(const Gmpq &a, const Gmpq &b)
{ return b < a; }

inline
bool
operator>=(const Gmpq &a, const Gmpq &b)
{ return ! (a < b); }

inline
bool
operator!=(const Gmpq &a, const Gmpq &b)
{ return ! (a == b); }


// mixed operators.
inline
bool
operator<(const Gmpq &a, int b)
{ return mpq_cmp_si(a.mpq(), b, 1) < 0; }

inline
bool
operator<(int a, const Gmpq &b)
{ return mpq_cmp_si(b.mpq(), a, 1) > 0; }

inline
bool
operator==(const Gmpq &a, int b)
{ return mpq_cmp_si(a.mpq(), b, 1) == 0; }

inline
bool
operator==(int a, const Gmpq &b)
{ return b == a; }

inline
bool
operator<=(const Gmpq &a, int b)
{ return ! (b < a); }

inline
bool
operator<=(int a, const Gmpq &b)
{ return ! (b < a); }

inline
bool
operator>(const Gmpq &a, int b)
{ return b < a; }

inline
bool
operator>(int a, const Gmpq &b)
{ return b < a; }

inline
bool
operator>=(const Gmpq &a, int b)
{ return ! (a < b); }

inline
bool
operator>=(int a, const Gmpq &b)
{ return ! (a < b); }

inline
bool
operator!=(const Gmpq &a, int b)
{ return ! (a == b); }

inline
bool
operator!=(int a, const Gmpq &b)
{ return ! (a == b); }



inline
Gmpq
Gmpq::operator-() const
{
    Gmpq Res;
    mpq_neg(Res.mpq(), mpq());
    return Res;
}

inline
Gmpq
operator+(const Gmpq &a, const Gmpq &b)
{
    Gmpq Res;
    mpq_add(Res.mpq(), a.mpq(), b.mpq());
    return Res;
}


inline
Gmpq&
Gmpq::operator+=(const Gmpq &z)
{
    *this = *this + z;
    return *this;
}


inline
Gmpq
operator-(const Gmpq &a, const Gmpq &b)
{
    Gmpq Res;
    mpq_sub(Res.mpq(), a.mpq(), b.mpq());
    return Res;
}


inline
Gmpq&
Gmpq::operator-=(const Gmpq &z)
{
    *this = *this - z;
    return *this;
}

inline
Gmpq
operator*(const Gmpq &a, const Gmpq &b)
{
    Gmpq Res;
    mpq_mul(Res.mpq(), a.mpq(), b.mpq());
    return Res;
}

inline
Gmpq&
Gmpq::operator*=(const Gmpq &z)
{
    *this = *this * z;
    return *this;
}


inline
Gmpq
operator/(const Gmpq &a, const Gmpq &b)
{
    CGAL_precondition(b != 0);
    Gmpq Res;
    mpq_div(Res.mpq(), a.mpq(), b.mpq());
    return Res;
}

inline
Gmpq&
Gmpq::operator/=(const Gmpq &z)
{
    *this = *this / z;
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
};

CGAL_END_NAMESPACE

#endif // CGAL_GMPQ_H
