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
// Author(s)     : Geert-Jan Giezeman
 

#ifndef CGAL_FLOAT_H
#define CGAL_FLOAT_H

#include <CGAL/basic.h>
#include <CGAL/tags.h>
#include <cmath>
#if defined(_MSC_VER) || defined(__BORLANDC__) || \
    defined(CGAL_MASK_FINITE_VALID) || defined(__PGI)
#  include <CGAL/IEEE_754_unions.h>
#endif
#ifdef __sgi
#  include <fp_class.h>
#endif

CGAL_BEGIN_NAMESPACE

template <> struct Number_type_traits<float> {
  typedef Tag_false Has_gcd;
  typedef Tag_true  Has_division;
  typedef Tag_false Has_sqrt;
};

inline
double
to_double(float f)
{ return static_cast<double>(f); }


inline 
std::pair<double,double>
to_interval(float f)
{
  return std::pair<double,double>(f, f);
}


#ifdef __sgi

inline
bool is_finite(float f)
{
    switch (fp_class_f(f)) {
    case FP_POS_NORM:
    case FP_NEG_NORM:
    case FP_POS_ZERO:
    case FP_NEG_ZERO:
    case FP_POS_DENORM:
    case FP_NEG_DENORM:
        return true;
    case FP_SNAN:
    case FP_QNAN:
    case FP_POS_INF:
    case FP_NEG_INF:
        return false;
    }
    return false; // NOT REACHED
}

inline
bool is_valid(float d)
{
    switch (fp_class_f(d)) {
    case FP_POS_NORM:
    case FP_NEG_NORM:
    case FP_POS_ZERO:
    case FP_NEG_ZERO:
    case FP_POS_INF:
    case FP_NEG_INF:
    case FP_POS_DENORM:
    case FP_NEG_DENORM:
        return true;
    case FP_SNAN:
    case FP_QNAN:
        return false;
    }
    return false; // NOT REACHED
}

#elif defined(_MSC_VER) || defined(__BORLANDC__) || \
      defined(CGAL_MASK_FINITE_VALID) || defined(__PGI)

#define CGAL_EXPONENT_FLOAT_MASK   0x7f800000
#define CGAL_MANTISSA_FLOAT_MASK   0x007fffff

inline
bool
is_finite_by_mask_float(unsigned int u)
{
  unsigned int e = u & CGAL_EXPONENT_FLOAT_MASK;
  return ( (e ^ CGAL_EXPONENT_FLOAT_MASK) != 0);
}

inline
bool
is_nan_by_mask_float(unsigned int u)
{
  if ( is_finite_by_mask_float(u) ) return false;
  // unsigned int m = u & CGAL_MANTISSA_FLOAT_MASK;
  return ( (u & CGAL_MANTISSA_FLOAT_MASK) != 0);
}

inline
bool
is_finite( const float& flt)
{
  float f = flt;
  IEEE_754_float* p = reinterpret_cast<IEEE_754_float*>(&f);
  return is_finite_by_mask_float( p->c );
}

inline
bool
is_valid( const float& flt)
{
  float f = flt;
  IEEE_754_float* p = reinterpret_cast<IEEE_754_float*>(&f);
  return !is_nan_by_mask_float( p->c );
}

#else

inline
bool
is_valid(float d)
{ return (d == d); }

inline
bool
is_finite(float d)
{ return (d == d) && (is_valid(d-d)); }

#endif

inline
io_Operator
io_tag(float)
{ return io_Operator(); }

CGAL_END_NAMESPACE

#endif // CGAL_FLOAT_H
