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
// Author(s)     : Stefan Schirra
 
#ifndef CGAL_LEDA_REAL_H
#define CGAL_LEDA_REAL_H

#include <CGAL/basic.h>
#include <CGAL/Number_type_traits.h>
#include <CGAL/Interval_arithmetic.h>

#include <utility>

#include <CGAL/LEDA_basic.h>
#include <LEDA/real.h>

CGAL_BEGIN_NAMESPACE

template <> struct Number_type_traits<leda_real> {
  typedef Tag_false Has_gcd;
  typedef Tag_true  Has_division;
  typedef Tag_true  Has_sqrt;
};

#ifndef CGAL_NO_NAMESPACE
inline
double
to_double(const leda_real & r)
{ return r.to_double(); }
#endif // CGAL_NO_NAMESPACE

inline
leda_real
sqrt(const leda_real & r)
{ return CGAL_LEDA_SCOPE::sqrt(r); }

inline
bool
is_finite(const leda_real &)
{ return true; }

inline
bool
is_valid(const leda_real &)
{ return true; }

inline
io_Operator
io_tag(const leda_real &)
{ return io_Operator(); }

#ifndef CGAL_CFG_NO_NAMESPACE
#ifndef CGAL_USE_ADL_FOR_NT
inline
Sign
sign(const leda_real& r)
{ return (Sign)CGAL_LEDA_SCOPE::sign(r); }

inline
Comparison_result
compare(const leda_real& r1, const leda_real& r2)
{
  int c = CGAL_LEDA_SCOPE::compare(r1,r2);
  return (c < 0) ? SMALLER : ((0 < c) ?  LARGER : EQUAL);
}
#endif // CGAL_USE_ADL_FOR_NT
#endif // CGAL_CFG_NO_NAMESPACE

inline
std::pair<double,double>
to_interval (const leda_real & z)
{
  Protect_FPU_rounding<true> P (CGAL_FE_TONEAREST);
  double approx = z.to_double();
  double rel_error = z.get_double_error();
  FPU_set_cw(CGAL_FE_UPWARD);
  Interval_nt_advanced ina(-rel_error,rel_error);
  ina += 1;
  ina *= approx;
  return ina.pair();
}

CGAL_END_NAMESPACE

#endif // CGAL_LEDA_REAL_H
