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
 

#ifndef CGAL_LEDA_BIGFLOAT_H
#define CGAL_LEDA_BIGFLOAT_H

#include <CGAL/basic.h>
#include <CGAL/Number_type_traits.h>
#include <CGAL/Interval_arithmetic.h>

#include <utility>

#include <CGAL/LEDA_basic.h>
#include <LEDA/bigfloat.h>

CGAL_BEGIN_NAMESPACE

template <> struct Number_type_traits<leda_bigfloat> {
  typedef Tag_false Has_gcd;
  typedef Tag_true  Has_division;
  typedef Tag_true  Has_sqrt;

  typedef Tag_true  Has_exact_ring_operations;
  typedef Tag_false Has_exact_division;
  typedef Tag_false Has_exact_sqrt;
};

#ifndef CGAL_CFG_NO_NAMESPACE
inline
double
to_double(const leda_bigfloat & b)
{ return CGAL_LEDA_SCOPE::to_double(b); }
#endif // CGAL_CFG_NO_NAMESPACE

inline
bool
is_finite(const leda_bigfloat & b)
{ return !( CGAL_LEDA_SCOPE::isInf(b) || CGAL_LEDA_SCOPE::isNaN(b) ); }

inline
bool
is_valid(const leda_bigfloat & b)
{ return !( CGAL_LEDA_SCOPE::isNaN(b) ); }

inline
io_Operator
io_tag(const leda_bigfloat &)
{ return io_Operator(); }

inline
std::pair<double,double>
to_interval (const leda_bigfloat & z)
{
  // assuming leda_bigfloat guarantee 1 bit error max
  Protect_FPU_rounding<true> P (CGAL_FE_TONEAREST);
  Interval_nt_advanced approx (CGAL_LEDA_SCOPE::to_double(z));
  FPU_set_cw(CGAL_FE_UPWARD);
  approx += Interval_nt<false>::smallest();
  return approx.pair();
}

CGAL_END_NAMESPACE

#endif // CGAL_LEDA_BIGFLOAT_H
