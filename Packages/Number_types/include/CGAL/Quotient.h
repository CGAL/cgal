// Copyright (c) 1999-2004  Utrecht University (The Netherlands),
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
// Author(s)     : Stefan Schirra, Sylvain Pion

// The template class Quotient<NT> is based on the LEDA class
// leda_rational written by Stefan Naeher and Christian Uhrig.
// It is basically a templated version with restricted functionality
// of the version of rational in LEDA release 3.3.
// The modification was done by Stefan.Schirra@mpi-sb.mpg.de

#ifndef CGAL_QUOTIENT_H
#define CGAL_QUOTIENT_H

#include <CGAL/basic.h>
#include <CGAL/MP_Float.h>
#include <utility>

#ifndef CGAL_CFG_NO_LOCALE
#  include <locale>
#else
#  include <cctype>
#endif

#include <CGAL/Interval_arithmetic.h>

CGAL_BEGIN_NAMESPACE

namespace CGALi {

  // Mini helper to prevent clashes when NT == int.
  template < typename T >
  struct Int_if_not_int { typedef int type; };

  template <>
  struct Int_if_not_int<int> { struct type{}; };

} // namespace CGALi

#define CGAL_int(T) typename CGALi::Int_if_not_int<T>::type

// Simplify the quotient numerator/denominator.
// Currently the default template doesn't do anything.
// This function is not documented as a number type requirement for now.
template < typename NT >
inline void
simplify_quotient(NT &, NT &) {}

template <class NT_>
class Quotient
{
 public:
  typedef NT_        NT;
  typedef Tag_false  Has_gcd;
  typedef Tag_true   Has_division;
  typedef typename Number_type_traits<NT_>::Has_sqrt  Has_sqrt;

  Quotient() : num(0), den(1) {}

  template <class T>
  Quotient(const T& n) : num(n), den(1) {}

  template <class T1, class T2>
  Quotient(const T1& n, const T2& d) : num(n), den(d)
  { CGAL_precondition( d != 0 ); }

  Quotient<NT>& operator+= (const Quotient<NT>& r);
  Quotient<NT>& operator-= (const Quotient<NT>& r);
  Quotient<NT>& operator*= (const Quotient<NT>& r);
  Quotient<NT>& operator/= (const Quotient<NT>& r);
  Quotient<NT>& operator+= (const NT& r);
  Quotient<NT>& operator-= (const NT& r);
  Quotient<NT>& operator*= (const NT& r);
  Quotient<NT>& operator/= (const NT& r);
  Quotient<NT>& operator+= (const CGAL_int(NT)& r);
  Quotient<NT>& operator-= (const CGAL_int(NT)& r);
  Quotient<NT>& operator*= (const CGAL_int(NT)& r);
  Quotient<NT>& operator/= (const CGAL_int(NT)& r);

  Quotient<NT>&    normalize();

  const NT&   numerator()   const { return num; }
  const NT&   denominator() const { return den; }

  void swap(Quotient &q)
  {
    using std::swap;
    swap(num, q.num);
    swap(den, q.den);
  }

 public:
  NT   num;
  NT   den;
};

template <class NT>
inline
void swap(Quotient<NT> &p, Quotient<NT> &q)
{
  p.swap(q);
}

#ifndef CGAL_NEW_NT_TRAITS
template <class NT>
Quotient<NT>
sqrt(const Quotient<NT> &q)
{
    CGAL_precondition(q > 0);
    return Quotient<NT>(CGAL_NTS sqrt(q.numerator()*q.denominator()),
	                q.denominator());
}
#endif

template <class NT>
CGAL_KERNEL_MEDIUM_INLINE
Quotient<NT>&
Quotient<NT>::normalize()
{
  if (num == den)
  {
      num = den = 1;
      return *this;
  }
  if (-num == den)
  {
      num = -1;
      den = 1;
      return *this;
  }
  NT ggt = gcd(num, den);
  if (ggt != 1 )
  {
      num /= ggt;
      den /= ggt;
  }
  return *this;
}

template <class NT>
CGAL_KERNEL_INLINE
Quotient<NT>&
Quotient<NT>::operator+= (const Quotient<NT>& r)
{
    num = num * r.den + r.num * den;
    den *= r.den;
    simplify_quotient(num, den);
    return *this;
}

template <class NT>
CGAL_KERNEL_INLINE
Quotient<NT>&
Quotient<NT>::operator-= (const Quotient<NT>& r)
{
    num = num * r.den - r.num * den;
    den *= r.den;
    simplify_quotient(num, den);
    return *this;
}

template <class NT>
CGAL_KERNEL_INLINE
Quotient<NT>&
Quotient<NT>::operator*= (const Quotient<NT>& r)
{
    num *= r.num;
    den *= r.den;
    simplify_quotient(num, den);
    return *this;
}

template <class NT>
CGAL_KERNEL_INLINE
Quotient<NT>&
Quotient<NT>::operator/= (const Quotient<NT>& r)
{
    CGAL_precondition( r.num != 0 );
    num *= r.den;
    den *= r.num;
    simplify_quotient(num, den);
    return *this;
}

template <class NT>
CGAL_KERNEL_INLINE
Quotient<NT>&
Quotient<NT>::operator+= (const NT& r)
{
    num += r * den;
    return *this;
}

template <class NT>
CGAL_KERNEL_INLINE
Quotient<NT>&
Quotient<NT>::operator-= (const NT& r)
{
    num -= r * den;
    return *this;
}

template <class NT>
CGAL_KERNEL_INLINE
Quotient<NT>&
Quotient<NT>::operator*= (const NT& r)
{
    num *= r;
    return *this;
}

template <class NT>
CGAL_KERNEL_INLINE
Quotient<NT>&
Quotient<NT>::operator/= (const NT& r)
{
    CGAL_precondition( r != 0 );
    den *= r;
    return *this;
}

template <class NT>
CGAL_KERNEL_INLINE
Quotient<NT>&
Quotient<NT>::operator+= (const CGAL_int(NT)& r)
{
    num += r * den;
    return *this;
}

template <class NT>
CGAL_KERNEL_INLINE
Quotient<NT>&
Quotient<NT>::operator-= (const CGAL_int(NT)& r)
{
    num -= r * den;
    return *this;
}

template <class NT>
CGAL_KERNEL_INLINE
Quotient<NT>&
Quotient<NT>::operator*= (const CGAL_int(NT)& r)
{
    num *= r;
    return *this;
}

template <class NT>
CGAL_KERNEL_INLINE
Quotient<NT>&
Quotient<NT>::operator/= (const CGAL_int(NT)& r)
{
    CGAL_precondition( r != 0 );
    den *= r;
    return *this;
}

#ifndef CGAL_NEW_NT_TRAITS
template <class NT>
CGAL_KERNEL_MEDIUM_INLINE
Comparison_result
quotient_cmp(const Quotient<NT>& x, const Quotient<NT>& y)
{
    // No assumptions on the sign of  den  are made

    // code assumes that SMALLER == - 1;
    CGAL_precondition( SMALLER == static_cast<Comparison_result>(-1) );

    int xsign = CGAL_NTS sign(x.num) * CGAL_NTS sign(x.den) ;
    int ysign = CGAL_NTS sign(y.num) * CGAL_NTS sign(y.den) ;
    if (xsign == 0) return static_cast<Comparison_result>(-ysign);
    if (ysign == 0) return static_cast<Comparison_result>(xsign);
    // now (x != 0) && (y != 0)
    int diff = xsign - ysign;
    if (diff == 0)
    {
        int msign = CGAL_NTS sign(x.den) * CGAL_NTS sign(y.den);
        NT leftop  = x.num * y.den * msign;
        NT rightop = y.num * x.den * msign;
        return CGAL_NTS compare(leftop, rightop);
    }
    else
    {
        return (xsign < ysign) ? SMALLER : LARGER;
    }
}

template <class NT>
inline
Comparison_result
compare(const Quotient<NT>& x, const Quotient<NT>& y)
{ return quotient_cmp(x, y); }
#endif

template <class NT>
std::ostream&
operator<<(std::ostream& s, const Quotient<NT>& r)
{
   return s << r.numerator() << "/" << r.denominator();
}

template <class NT>
std::istream&
operator>>(std::istream& in, Quotient<NT>& r)
{
  /* format  num/den  or simply  num  */

  char c = 0;

#ifndef CGAL_CFG_NO_LOCALE
  while (in.get(c) && std::isspace(c, std::locale::classic() ));
#else
  while (in.get(c) && CGAL_CLIB_STD::isspace(c));
#endif // CGAL_CFG_NO_LOCALE
  if ( !in ) return in;
  in.putback(c);

  NT num;
  NT den(1);
  in >> num;

#ifndef CGAL_CFG_NO_LOCALE
  while (in.get(c) && std::isspace(c, std::locale::classic() ));
#else
  while (in.get(c) && CGAL_CLIB_STD::isspace(c));
#endif // CGAL_CFG_NO_LOCALE
  if (( in ) && ( c == '/'))
  {
#ifndef CGAL_CFG_NO_LOCALE
      while (in.get(c) && std::isspace(c, std::locale::classic() ));
#else
      while (in.get(c) && CGAL_CLIB_STD::isspace(c));
#endif // CGAL_CFG_NO_LOCALE
      CGAL_assertion( in );
      in.putback(c);
      in >> den;
  }
  else
  {
      in.putback(c);
      if ( in.eof() ) in.clear();
  }
  r = Quotient<NT>( num, den);
  return in;
}

#ifndef CGAL_NEW_NT_TRAITS
template <class NT>
inline
io_Operator
io_tag(const Quotient<NT>&)
{ return io_Operator(); }
#endif

template <class NT>
CGAL_KERNEL_INLINE
Quotient<NT>
operator+(const Quotient<NT>& x, const Quotient<NT>& y)
{
  Quotient<NT> z = x;
  return z += y;
}

template <class NT>
CGAL_KERNEL_INLINE
Quotient<NT>
operator-(const Quotient<NT>& x, const Quotient<NT>& y)
{ return (Quotient<NT>(x) -= y); }

template <class NT>
CGAL_KERNEL_INLINE
Quotient<NT>
operator*(const Quotient<NT>& x, const Quotient<NT>& y)
{
  Quotient<NT> z = x;
  return z *= y;
}

template <class NT>
CGAL_KERNEL_INLINE
Quotient<NT>
operator/(const Quotient<NT>& x, const Quotient<NT>& y)
{
  Quotient<NT> z = x;
  return z /= y;
}

template <class NT>
inline
Quotient<NT>
operator-(const Quotient<NT>& x)
{ return Quotient<NT>(-x.num,x.den); }

template <class NT>
CGAL_KERNEL_INLINE
Quotient<NT>
operator+(const NT& x, const Quotient<NT>& y)
{
  Quotient<NT> z(x);
  return z += y;
}

template <class NT>
CGAL_KERNEL_INLINE
Quotient<NT>
operator+(const Quotient<NT>& x, const NT& y)
{
  Quotient<NT> z = x;
  return z += y;
}

template <class NT>
CGAL_KERNEL_INLINE
Quotient<NT>
operator+(const Quotient<NT>& x, const CGAL_int(NT)& y)
{
  Quotient<NT> z = x;
  return z += NT(y);
}

template <class NT>
CGAL_KERNEL_INLINE
Quotient<NT>
operator+(const CGAL_int(NT)& x, const Quotient<NT>& y)
{ return y + x; }

template <class NT>
CGAL_KERNEL_INLINE
Quotient<NT>
operator-(const NT& x, const Quotient<NT>& y)
{
  Quotient<NT> z(x);
  return z -= y;
}

template <class NT>
CGAL_KERNEL_INLINE
Quotient<NT>
operator-(const Quotient<NT>& x, const NT& y)
{
  Quotient<NT> z = x;
  return z -= y;
}

template <class NT>
CGAL_KERNEL_INLINE
Quotient<NT>
operator-(const Quotient<NT>& x, const CGAL_int(NT)& y)
{
  Quotient<NT> z = x;
  return z -= NT(y);
}

template <class NT>
CGAL_KERNEL_INLINE
Quotient<NT>
operator-(const CGAL_int(NT)& x, const Quotient<NT>& y)
{
  Quotient<NT> z = x;
  return z -= y;
}

template <class NT>
CGAL_KERNEL_INLINE
Quotient<NT>
operator*(const NT& x, const Quotient<NT>& y)
{
  Quotient<NT> z(x);
  return z *= y;
}

template <class NT>
CGAL_KERNEL_INLINE
Quotient<NT>
operator*(const Quotient<NT>& x, const NT& y)
{
  Quotient<NT> z = x;
  return z *= y;
}

template <class NT>
CGAL_KERNEL_INLINE
Quotient<NT>
operator*(const Quotient<NT>& x, const CGAL_int(NT)& y)
{
  Quotient<NT> z = x;
  return z *= NT(y);
}

template <class NT>
inline
Quotient<NT>
operator*(const CGAL_int(NT)& x, const Quotient<NT>& y)
{ return y*x; }

template <class NT>
CGAL_KERNEL_INLINE
Quotient<NT>
operator/(const NT& x, const Quotient<NT>& y)
{
  Quotient<NT> z(x) ;
  return z /= y;
}

template <class NT>
CGAL_KERNEL_INLINE
Quotient<NT>
operator/(const Quotient<NT>& x, const NT& y)
{
  Quotient<NT> z = x;
  return z /= y;
}

template <class NT>
CGAL_KERNEL_INLINE
Quotient<NT>
operator/(const Quotient<NT>& x, const CGAL_int(NT)& y)
{
  Quotient<NT> z = x;
  return z /= NT(y);
}

template <class NT>
CGAL_KERNEL_INLINE
Quotient<NT>
operator/(const CGAL_int(NT)& x, const Quotient<NT>& y)
{
  Quotient<NT> z = x;
  return z /= y;
}

template <class NT>
CGAL_KERNEL_INLINE
NT
quotient_truncation(const Quotient<NT>& r)
{ return (r.num / r.den); }



template <class NT>
CGAL_KERNEL_INLINE
bool
operator==(const Quotient<NT>& x, const Quotient<NT>& y)
{ return x.num * y.den == x.den * y.num; }

template <class NT>
CGAL_KERNEL_INLINE
bool
operator==(const Quotient<NT>& x, const NT& y)
{ return x.den * y == x.num; }

template <class NT>
inline
bool
operator==(const NT& x, const Quotient<NT>& y)
{ return y == x; }

template <class NT>
CGAL_KERNEL_INLINE
bool
operator==(const CGAL_int(NT) & x, const Quotient<NT>& y)
{ return y.den * x == y.num; }

template <class NT>
inline
bool
operator==(const Quotient<NT>& x, const CGAL_int(NT) & y)
{ return y == x; }


template <class NT>
inline
bool
operator!=(const Quotient<NT>& x, const Quotient<NT>& y)
{ return ! (x == y); }

template <class NT>
inline
bool
operator!=(const Quotient<NT>& x, const NT& y)
{ return ! (x == y); }

template <class NT>
inline
bool
operator!=(const NT& x, const Quotient<NT>& y)
{ return ! (x == y); }

template <class NT>
inline
bool
operator!=(const CGAL_int(NT) & x, const Quotient<NT>& y)
{ return ! (x == y); }

template <class NT>
inline
bool
operator!=(const Quotient<NT>& x, const CGAL_int(NT) & y)
{ return ! (x == y); }


template <class NT>
CGAL_KERNEL_INLINE
bool
operator<(const Quotient<NT>& x, const Quotient<NT>& y)
{
  return quotient_cmp(x,y) == SMALLER;
}

template <class NT>
CGAL_KERNEL_INLINE
bool
operator<(const Quotient<NT>& x, const NT& y)
{
  return quotient_cmp(x,Quotient<NT>(y)) == SMALLER;
}

template <class NT>
CGAL_KERNEL_INLINE
bool
operator<(const NT& x, const Quotient<NT>& y)
{
  return quotient_cmp(Quotient<NT>(x),y) == SMALLER;
}

template <class NT>
CGAL_KERNEL_INLINE
bool
operator<(const CGAL_int(NT)& x, const Quotient<NT>& y)
{
  return quotient_cmp(Quotient<NT>(x),y) == SMALLER;
}

template <class NT>
CGAL_KERNEL_INLINE
bool
operator<(const Quotient<NT>& x, const CGAL_int(NT)& y)
{
  return quotient_cmp(x,Quotient<NT>(y)) == SMALLER;
}


template <class NT>
inline
bool
operator>(const Quotient<NT>& x, const Quotient<NT>& y)
{ return y < x; }

template <class NT>
inline
bool
operator>(const Quotient<NT>& x, const NT& y)
{ return y < x; }

template <class NT>
inline
bool
operator>(const NT& x, const Quotient<NT>& y)
{ return y < x; }

template <class NT>
inline
bool
operator>(const Quotient<NT>& x, const CGAL_int(NT)& y)
{ return y < x; }

template <class NT>
inline
bool
operator>(const CGAL_int(NT)& x, const Quotient<NT>& y)
{ return y < x; }



template <class NT>
inline
bool
operator<=(const Quotient<NT>& x, const Quotient<NT>& y)
{ return ! (y < x); }

template <class NT>
inline
bool
operator<=(const Quotient<NT>& x, const NT& y)
{ return ! (y < x); }

template <class NT>
inline
bool
operator<=(const NT& x, const Quotient<NT>& y)
{ return ! (y < x); }

template <class NT>
inline
bool
operator<=(const Quotient<NT>& x, const CGAL_int(NT)& y)
{ return ! (y < x); }

template <class NT>
inline
bool
operator<=(const CGAL_int(NT)& x, const Quotient<NT>& y)
{ return ! (y < x); }



template <class NT>
inline
bool
operator>=(const Quotient<NT>& x, const Quotient<NT>& y)
{ return ! (x < y); }

template <class NT>
inline
bool
operator>=(const Quotient<NT>& x, const NT& y)
{ return ! (x < y); }

template <class NT>
inline
bool
operator>=(const NT& x, const Quotient<NT>& y)
{ return ! (x < y); }

template <class NT>
inline
bool
operator>=(const Quotient<NT>& x, const CGAL_int(NT)& y)
{ return ! (x < y); }

template <class NT>
inline
bool
operator>=(const CGAL_int(NT)& x, const Quotient<NT>& y)
{ return ! (x < y); }


#ifdef CGAL_NEW_NT_TRAITS

template<class NT> struct Number_type_traits_base;

template<class NT>
class Number_type_traits_base< Quotient<NT> >
  : public CGALi::Default_field_number_type_traits< Quotient<NT> >
{
private:
  typedef Quotient<NT>  QNT;

public:
  typedef typename Number_type_traits<NT>::Has_exact_ring_operations
  Has_exact_ring_operations;

  typedef Has_exact_ring_operations   Has_exact_division;

  typedef Tag_true Has_sqrt;
  typedef typename Number_type_traits<NT>::Has_exact_sqrt
  Has_exact_sqrt;

  typedef Tag_true Has_rational_traits;

  static inline double to_double(const QNT& q) {
    if (q.num == 0 )
      return 0;

    double nd = CGAL::to_double( q.num );

    if (q.den == 1 )
      return nd;

    double dd = CGAL::to_double( q.den );

    if ( is_finite( q.den ) && is_finite( q.num ) )
      return nd/dd;

    if ( CGAL::abs(q.num) > CGAL::abs(q.den) )
      {
	NT  nt_div = q.num / q.den;
	double divd = CGAL::to_double(nt_div);
	if ( divd >= CGAL_CLIB_STD::ldexp(1.0,53) )
	  { return divd; }
      }
    if ( CGAL::abs(q.num) < CGAL::abs(q.den) )
      { return 1.0 / CGAL::to_double( NT(1) / q ); }

    return nd/dd;
  }

  static inline std::pair<double,double>
  to_interval (const QNT& z) {
    Interval_nt<> quot = Interval_nt<>(CGAL::to_interval(z.numerator())) /
		         Interval_nt<>(CGAL::to_interval(z.denominator()));
    return std::make_pair(quot.inf(), quot.sup());
  }

  static inline bool is_valid(const QNT& q) {
    return is_valid(q.num) && is_valid(q.den);
  }

  static inline bool is_finite(const QNT& q) {
    return is_finite(q.num) && is_finite(q.den);
  }

  static inline QNT sqrt(const QNT &q) {
    CGAL_precondition(q > 0);
    return QNT(CGAL::sqrt(q.numerator()*q.denominator()),
	       q.denominator());
  }

private:
  static inline Comparison_result
  quotient_cmp(const QNT& x, const QNT& y)
  {
    // No assumptions on the sign of  den  are made

    // code assumes that SMALLER == - 1;
    CGAL_precondition( SMALLER == static_cast<Comparison_result>(-1) );

    int xsign = CGAL::sign(x.num) * CGAL::sign(x.den) ;
    int ysign = CGAL::sign(y.num) * CGAL::sign(y.den) ;
    if (xsign == 0) return static_cast<Comparison_result>(-ysign);
    if (ysign == 0) return static_cast<Comparison_result>(xsign);
    // now (x != 0) && (y != 0)
    int diff = xsign - ysign;
    if (diff == 0)
      {
	int msign = CGAL::sign(x.den) * CGAL::sign(y.den);
        NT leftop  = x.num * y.den * msign;
        NT rightop = y.num * x.den * msign;
        return CGAL::compare(leftop, rightop);
      }
    else
      {
        return (xsign < ysign) ? SMALLER : LARGER;
      }
  }

public:
  static inline Comparison_result
  compare(const QNT& x, const QNT& y) {
    return quotient_cmp(x, y);
  }

};

template<class NT>
struct Number_type_traits< Quotient<NT> >
  : public Number_type_traits_base< Quotient<NT> > {};


template<>
struct Number_type_traits< Quotient<MP_Float> >
  : public Number_type_traits_base< Quotient<MP_Float> >
{
  using std::pair;

  static double to_double(const Quotient<MP_Float> &q)
  {
    pair<double, int> n = to_double_exp(q.numerator());
    pair<double, int> d = to_double_exp(q.denominator());
    double scale = CGAL_CLIB_STD::ldexp(1.0, n.second - d.second);
    return (n.first / d.first) * scale;
  }

  // FIXME : This function deserves proper testing...
  static pair<double,double>
  to_interval(const Quotient<MP_Float> &q)
  {
    pair<pair<double, double>, int> n = to_interval_exp(q.numerator());
    pair<pair<double, double>, int> d = to_interval_exp(q.denominator());
    double scale = CGAL_CLIB_STD::ldexp(1.0, n.second - d.second);
    return ((Interval_nt<>(n.first) / Interval_nt<>(d.first)) * scale).pair();
  }
};

#else // CGAL_NEW_NT_TRAITS

template <class NT>
double
to_double(const Quotient<NT>& q)   /* TODO */
{
  if (q.num == 0 )
    return 0;

  double nd = CGAL_NTS to_double( q.num );

  if (q.den == 1 )
    return nd;

  double dd = CGAL_NTS to_double( q.den );

  if ( is_finite( q.den ) && is_finite( q.num ) )
    return nd/dd;

  if ( CGAL_NTS abs(q.num) > CGAL_NTS abs(q.den) )
  {
      NT  nt_div = q.num / q.den;
      double divd = CGAL_NTS to_double(nt_div);
      if ( divd >= CGAL_CLIB_STD::ldexp(1.0,53) )
      { return divd; }
  }
  if ( CGAL_NTS abs(q.num) < CGAL_NTS abs(q.den) )
  { return 1.0 / CGAL_NTS to_double( NT(1) / q ); }

  return nd/dd;
}

template <class RT>
std::pair<double,double>
to_interval (const Quotient<RT>& z)
{
    Interval_nt<> quot = Interval_nt<>(CGAL_NTS to_interval(z.numerator())) /
		         Interval_nt<>(CGAL_NTS to_interval(z.denominator()));
    return std::make_pair(quot.inf(), quot.sup());
}

template <class NT>
bool
is_valid(const Quotient<NT>& q)
{ return is_valid(q.num) && is_valid(q.den); }

template <class NT>
bool
is_finite(const Quotient<NT>& q)
{ return is_finite(q.num) && is_finite(q.den); }

#endif // CGAL_NEW_NT_TRAITS

template <class NT>
inline
const NT&
denominator(const Quotient<NT>& q)
{ return q.den ; }

template <class NT>
inline
const NT&
numerator(const Quotient<NT>& q)
{ return q.num ; }

/*
template <class NT>
NT
gcd(const NT&, const NT&)
{ return NT(1); }
*/

#undef CGAL_int

CGAL_END_NAMESPACE

#endif  // CGAL_QUOTIENT_H
