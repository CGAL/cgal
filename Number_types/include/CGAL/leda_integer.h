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
 
#ifndef CGAL_LEDA_INTEGER_H
#define CGAL_LEDA_INTEGER_H

#include <CGAL/basic.h>
#include <CGAL/Number_type_traits.h>
#include <CGAL/Interval_nt.h>

#include <utility>

#include <CGAL/LEDA_basic.h>
#include <LEDA/integer.h>

CGAL_BEGIN_NAMESPACE

template <> struct Number_type_traits<leda_integer> {
  typedef Tag_true  Has_gcd;
  typedef Tag_true  Has_division;
  typedef Tag_true  Has_sqrt;

  typedef Tag_true  Has_exact_ring_operations;
  typedef Tag_false Has_exact_division;
  typedef Tag_false Has_exact_sqrt;
};

inline
double
to_double(const leda_integer & i)
{ return i.to_double(); }

inline
leda_integer
sqrt(const leda_integer & i)
{ return CGAL_LEDA_SCOPE::sqrt(i); }

inline
bool
is_finite(const leda_integer &)
{ return true; }

inline
bool
is_valid(const leda_integer &)
{ return true; }

inline
io_Operator
io_tag(const leda_integer &)
{ return io_Operator(); }

inline
Sign
sign(const leda_integer& n)
{ return (Sign) CGAL_LEDA_SCOPE::sign(n); }

inline
leda_integer
div( const leda_integer& n1, const leda_integer& n2)
{ 
  return n1 / n2;
}

// missing mixed operators
inline
bool
operator==(int a, const leda_integer& b)
{ return b == a; }

inline
bool
operator!=(int a, const leda_integer& b)
{ return b != a; }


inline
std::pair<double,double>
to_interval (const leda_integer & n)
{
  Protect_FPU_rounding<true> P (CGAL_FE_TONEAREST);
  double cn = CGAL::to_double(n);
  leda_integer pn = ( n>0 ? n : -n);
  if ( pn.iszero() || log(pn) < 53 )
      return to_interval(cn);
  else {
    FPU_set_cw(CGAL_FE_UPWARD);
    Interval_nt_advanced ina(cn);
    ina += Interval_nt_advanced::smallest();
    return ina.pair();
  }
}

inline
leda_integer
gcd( const leda_integer& n1, const leda_integer& n2)
{ 
  return CGAL_LEDA_SCOPE::gcd(n1, n2);
}

CGAL_END_NAMESPACE

#endif // CGAL_LEDA_INTEGER_H
