// Copyright (c) 1999-2005  Utrecht University (The Netherlands),
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

// The include is done before the protect macro on purpose, because
// of a cyclic dependency.
#include <CGAL/basic.h>

#ifndef CGAL_QUOTIENT_H
#define CGAL_QUOTIENT_H

#include <utility>

#ifndef CGAL_CFG_NO_LOCALE
#  include <locale>
#else
#  include <cctype>
#endif

#include <CGAL/Interval_arithmetic.h>
#include <CGAL/Number_type_traits.h>
#include <CGAL/Kernel/mpl.h>
#include <CGAL/Binary_operator_result.h>
#include <CGAL/Quotient_fwd.h>

#include <boost/operators.hpp>

CGAL_BEGIN_NAMESPACE

#define CGAL_int(T)    typename First_if_different<int,    T>::Type
#define CGAL_double(T) typename First_if_different<double, T>::Type

// Simplify the quotient numerator/denominator.
// Currently the default template doesn't do anything.
// This function is not documented as a number type requirement for now.
template < typename NT >
inline void
simplify_quotient(NT &, NT &) {}

template <class NT_>
class Quotient
  : boost::ordered_field_operators1< Quotient<NT_>
  , boost::ordered_field_operators2< Quotient<NT_>, NT_
  , boost::ordered_field_operators2< Quotient<NT_>, CGAL_int(NT_)
    > > >
{
 public:
  typedef NT_        NT;
  typedef Tag_false  Has_gcd;
  typedef Tag_true   Has_division;
  typedef typename Number_type_traits<NT_>::Has_sqrt  Has_sqrt;

  typedef Tag_true   Has_exact_division;
  typedef typename Number_type_traits<NT_>::Has_exact_sqrt Has_exact_sqrt;
  typedef typename Number_type_traits<NT_>::Has_exact_ring_operations
  Has_exact_ring_operations;

  Quotient()
    : num(0), den(1) {}

  Quotient(const NT& n)
    : num(n), den(1) {}

  Quotient(const CGAL_double(NT) & n)
    : num(n), den(1) {}

  Quotient(const CGAL_int(NT) & n)
    : num(n), den(1) {}

  template <class T>
  explicit Quotient(const T& n) : num(n), den(1) {}

  template <class T>
  Quotient(const Quotient<T>& n) : num(n.numerator()), den(n.denominator()) {}

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


template < typename NT >
struct Binary_operator_result < Quotient<NT>, NT >
{ typedef Quotient<NT>  type; };

template < typename NT >
struct Binary_operator_result < NT, Quotient<NT> >
{ typedef Quotient<NT>  type; };


template <class NT>
inline
void swap(Quotient<NT> &p, Quotient<NT> &q)
{
  p.swap(q);
}

template <class NT>
Quotient<NT>
sqrt(const Quotient<NT> &q)
{
    CGAL_precondition(q > 0);
    return Quotient<NT>(CGAL_NTS sqrt(q.numerator()*q.denominator()),
	                q.denominator());
}

template <class NT>
CGAL_MEDIUM_INLINE
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
CGAL_MEDIUM_INLINE
Quotient<NT>&
Quotient<NT>::operator+= (const Quotient<NT>& r)
{
    num = num * r.den + r.num * den;
    den *= r.den;
    simplify_quotient(num, den);
    return *this;
}

template <class NT>
CGAL_MEDIUM_INLINE
Quotient<NT>&
Quotient<NT>::operator-= (const Quotient<NT>& r)
{
    num = num * r.den - r.num * den;
    den *= r.den;
    simplify_quotient(num, den);
    return *this;
}

template <class NT>
CGAL_MEDIUM_INLINE
Quotient<NT>&
Quotient<NT>::operator*= (const Quotient<NT>& r)
{
    num *= r.num;
    den *= r.den;
    simplify_quotient(num, den);
    return *this;
}

template <class NT>
CGAL_MEDIUM_INLINE
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
CGAL_MEDIUM_INLINE
Quotient<NT>&
Quotient<NT>::operator+= (const NT& r)
{
    num += r * den;
    return *this;
}

template <class NT>
CGAL_MEDIUM_INLINE
Quotient<NT>&
Quotient<NT>::operator-= (const NT& r)
{
    num -= r * den;
    return *this;
}

template <class NT>
CGAL_MEDIUM_INLINE
Quotient<NT>&
Quotient<NT>::operator*= (const NT& r)
{
    num *= r;
    return *this;
}

template <class NT>
CGAL_MEDIUM_INLINE
Quotient<NT>&
Quotient<NT>::operator/= (const NT& r)
{
    CGAL_precondition( r != 0 );
    den *= r;
    return *this;
}

template <class NT>
CGAL_MEDIUM_INLINE
Quotient<NT>&
Quotient<NT>::operator+= (const CGAL_int(NT)& r)
{
    num += r * den;
    return *this;
}

template <class NT>
CGAL_MEDIUM_INLINE
Quotient<NT>&
Quotient<NT>::operator-= (const CGAL_int(NT)& r)
{
    num -= r * den;
    return *this;
}

template <class NT>
CGAL_MEDIUM_INLINE
Quotient<NT>&
Quotient<NT>::operator*= (const CGAL_int(NT)& r)
{
    num *= r;
    return *this;
}

template <class NT>
CGAL_MEDIUM_INLINE
Quotient<NT>&
Quotient<NT>::operator/= (const CGAL_int(NT)& r)
{
    CGAL_precondition( r != 0 );
    den *= r;
    return *this;
}

template <class NT>
CGAL_MEDIUM_INLINE
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

template <class NT>
inline
io_Operator
io_tag(const Quotient<NT>&)
{ return io_Operator(); }


template <class NT>
inline
Quotient<NT>
operator-(const Quotient<NT>& x)
{ return Quotient<NT>(-x.num,x.den); }


template <class NT>
CGAL_MEDIUM_INLINE
NT
quotient_truncation(const Quotient<NT>& r)
{ return (r.num / r.den); }



template <class NT>
CGAL_MEDIUM_INLINE
bool
operator==(const Quotient<NT>& x, const Quotient<NT>& y)
{ return x.num * y.den == x.den * y.num; }

template <class NT>
CGAL_MEDIUM_INLINE
bool
operator==(const Quotient<NT>& x, const NT& y)
{ return x.den * y == x.num; }

template <class NT>
inline
bool
operator==(const Quotient<NT>& x, const CGAL_int(NT) & y)
{ return x.den * y == x.num; }



template <class NT>
CGAL_MEDIUM_INLINE
bool
operator<(const Quotient<NT>& x, const Quotient<NT>& y)
{
  return quotient_cmp(x,y) == SMALLER;
}

template <class NT>
CGAL_MEDIUM_INLINE
bool
operator<(const Quotient<NT>& x, const NT& y)
{
  return quotient_cmp(x,Quotient<NT>(y)) == SMALLER;
}

template <class NT>
CGAL_MEDIUM_INLINE
bool
operator<(const Quotient<NT>& x, const CGAL_int(NT)& y)
{
  return quotient_cmp(x,Quotient<NT>(y)) == SMALLER;
}


template <class NT>
inline
bool
operator>(const Quotient<NT>& x, const NT& y)
{ return quotient_cmp(x,Quotient<NT>(y)) == LARGER; }

template <class NT>
inline
bool
operator>(const Quotient<NT>& x, const CGAL_int(NT)& y)
{ return quotient_cmp(x, Quotient<NT>(y)) == LARGER; }


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

// The min/max are functions are needed since LEDA defines template
// min/max functions which clash with std::min/max with ADL.
template <class NT>
inline
const Quotient<NT>&
min(const Quotient<NT>& p, const Quotient<NT>& q)
{
  return std::min(p, q);
}
template <class NT>
inline
const Quotient<NT>&
max(const Quotient<NT>& p, const Quotient<NT>& q)
{
  return std::max(p, q);
}

/*
template <class NT>
NT
gcd(const NT&, const NT&)
{ return NT(1); }
*/

#undef CGAL_double
#undef CGAL_int

// Rational traits
template < class NT >
struct Rational_traits< Quotient<NT> >
{
  typedef NT RT;

  const RT & numerator   (const Quotient<NT>& r) const { return r.numerator(); }
  const RT & denominator (const Quotient<NT>& r) const { return r.denominator(); }
  
  Quotient<NT> make_rational(const RT & n, const RT & d) const
  { return Quotient<NT>(n, d); } 
};

CGAL_END_NAMESPACE

#endif  // CGAL_QUOTIENT_H
