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

#ifndef CGAL_NUMBER_UTILS_H
#define CGAL_NUMBER_UTILS_H

#include <CGAL/config.h>
#include <CGAL/enum.h>
#include <algorithm>

CGAL_BEGIN_NAMESPACE

template <class NT>
inline
bool
is_zero(const NT& x)
{ return x == 0; }

template <class NT>
inline
bool
is_one(const NT& x)
{ return x == 1; }

template <class NT>
inline
bool
is_negative(const NT& x)
{ return x < 0; }

template <class NT>
inline
bool
is_positive(const NT& x)
{ return 0 < x; }

template <class NT>
inline
Sign
sign(const NT& x)
{ return (x < 0) ? NEGATIVE : (0 < x) ? POSITIVE : ZERO; }

template <class NT>
inline
NT
abs(const NT& x)
{
  if (x < 0)
    return -x;
  return x;
}

template <class NT1, class NT2>
inline
Comparison_result
compare(const NT1& n1, const NT2& n2)
{ return (n1 < n2) ? SMALLER : (n2 < n1) ? LARGER : EQUAL; }

template <class NT>
inline
NT
square( const NT& n)
{ return n*n; }

template <class NT>
inline
NT
gcd( const NT& n1, const NT& n2)
{
  CGAL_precondition(!CGAL_NTS is_zero(n2));
  NT x = CGAL_NTS abs(n1);
  NT y = CGAL_NTS abs(n2);
  do {
    x %= y;
    if (CGAL_NTS is_zero(x)) return y;
    y %= x;
  } while (CGAL_NTS is_positive(y));
  return x;
}

// for min and max see <CGAL/number_type_basic.h>

CGAL_END_NAMESPACE

#endif // CGAL_NUMBER_UTILS_H
