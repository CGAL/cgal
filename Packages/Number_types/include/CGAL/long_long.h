// Copyright (c) 1999,2001  Utrecht University (The Netherlands),
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

// ISO C++ does not support `long long', but ISO C does, which means the next
// revision of ISO C++ probably will too.  However, currently, g++ -pedantic
// produces a warning so we don't include this file by default.

#ifndef CGAL_LONG_LONG_H
#define CGAL_LONG_LONG_H

#include <CGAL/basic.h>
#include <CGAL/Interval_arithmetic.h>

CGAL_BEGIN_NAMESPACE

template <> struct Number_type_traits<long long int> {
  typedef Tag_true   Has_gcd;
  typedef Tag_true   Has_division;
  typedef Tag_false  Has_sqrt;

  typedef Tag_false  Has_exact_ring_operations;
  typedef Tag_false  Has_exact_division;
  typedef Tag_false  Has_exact_sqrt;
};

inline
long long int
div(long long int i1, long long int i2)
{ return i1 / i2; }

inline
double
to_double(long long int i)
{ return (double)i; }

inline
bool
is_finite(long long int)
{ return true; }

inline
bool
is_valid(long long int)
{ return true; }

inline
std::pair<double,double>
to_interval (const long long & z)
{
  Protect_FPU_rounding<true> P(CGAL_FE_TONEAREST);
  Interval_nt<false> approx ((double) z);
  FPU_set_cw(CGAL_FE_UPWARD);
  approx += Interval_nt<false>::smallest();
  return approx.pair();
}


#if (defined(__sparc__) || defined(__sparc) || defined(sparc)) || \
    (defined(__sgi__)   || defined(__sgi)   || defined(sgi)) || \
    (defined(__i386__)  || defined(__i386)  || defined(i386)) || \
    (defined(__ppc__)   || defined(__ppc)   || defined(ppc)) || \
    (defined(__powerpc__) || defined(__powerpc) || defined(powerpc))
typedef  long long int           Integer64;
#define CGAL_HAS_INTEGER64
#endif

CGAL_END_NAMESPACE

#endif // CGAL_LONG_LONG_H
