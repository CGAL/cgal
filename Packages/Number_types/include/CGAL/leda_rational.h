// Copyright (c) 1999  Utrecht University (The Netherlands),
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
// Author(s)     : Andreas Fabri


#ifndef CGAL_LEDA_RATIONAL_H
#define CGAL_LEDA_RATIONAL_H

#include <CGAL/basic.h>
#include <CGAL/Number_type_traits.h>
#include <CGAL/Interval_nt.h>

#include <utility>

#include <CGAL/LEDA_basic.h>
#include <LEDA/rational.h>

CGAL_BEGIN_NAMESPACE

template <>
struct Number_type_traits<leda_rational> {
  typedef Tag_false Has_gcd;
  typedef Tag_true  Has_division;
  typedef Tag_false Has_sqrt;

  typedef Tag_true  Has_exact_ring_operations;
  typedef Tag_true  Has_exact_division;
  typedef Tag_false Has_exact_sqrt;
};

template <>
struct Rational_traits<leda_rational> {
  typedef leda_integer RT;
  RT numerator   (const leda_rational & r) const { return r.numerator(); }
  RT denominator (const leda_rational & r) const { return r.denominator(); }
  leda_rational make_rational(const RT & n, const RT & d) const
  { return leda_rational(n, d); }
  leda_rational make_rational(const leda_rational & n,
                              const leda_rational & d) const
  { return n / d; }
};

#ifndef CGAL_NO_NAMESPACE
inline
double
to_double(const leda_rational &r)
{ return r.to_double(); }
#endif // CGAL_NO_NAMESPACE

inline
bool
is_finite(const leda_rational &)
{ return true; }

inline
bool
is_valid(const leda_rational &)
{ return true; }

inline
io_Operator
io_tag(const leda_rational &)
{ return io_Operator(); }

inline
Sign
sign(const leda_rational& r)
{ return (Sign) CGAL_LEDA_SCOPE::sign(r); }

inline
std::pair<double,double>
to_interval (const leda_rational & z)
{
  // There's no guarantee about the error of to_double(), so I add 3 ulps...
  Protect_FPU_rounding<true> P (CGAL_FE_TONEAREST);
  Interval_nt_advanced approx (z.to_double());
  FPU_set_cw(CGAL_FE_UPWARD);

  approx += Interval_nt<false>::smallest();
  approx += Interval_nt<false>::smallest();
  approx += Interval_nt<false>::smallest();
  return approx.pair();
}

CGAL_END_NAMESPACE

#endif  // CGAL_LEDA_RATIONAL_H
