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

#ifndef CGAL_DOUBLE_H
#define CGAL_DOUBLE_H

#include <CGAL/basic.h>
#include <utility>
#include <cmath>
#ifdef CGAL_CFG_IEEE_754_BUG
#  include <CGAL/IEEE_754_unions.h>
#  define CGAL_EXPONENT_DOUBLE_MASK   0x7ff00000
#  define CGAL_MANTISSA_DOUBLE_MASK   0x000fffff
#endif
#ifdef __sgi
#  include <fp_class.h>
#endif

CGAL_BEGIN_NAMESPACE

#ifdef CGAL_NEW_NT_TRAITS

template<>
class Number_type_traits<double>
  : public CGALi::Default_field_number_type_traits<double>
{
public:
  typedef Tag_true   Has_sqrt;
  typedef Tag_false  Has_exact_division;
  typedef Tag_false  Has_exact_sqrt;

  typedef Tag_false  Has_rational_traits;
  typedef Tag_false  Has_simplify;

  static inline double to_double(double d) { return d; }

  static inline std::pair<double,double> to_interval(double d) {
    return std::make_pair(d, d);
  }

  static inline double sqrt(double d) { return std::sqrt(d); }

#ifdef __sgi
  static inline bool is_finite(double d) {
    switch (fp_class_d(d)) {
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

  static inline bool is_valid(double d) {
    switch (fp_class_d(d)) {
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

#elif defined CGAL_CFG_IEEE_754_BUG
private:
  static inline bool is_finite_by_mask_double(unsigned int h) {
    unsigned int e = h & CGAL_EXPONENT_DOUBLE_MASK;
    return ( ( e ^ CGAL_EXPONENT_DOUBLE_MASK ) != 0 );
  }

  static inline bool is_nan_by_mask_double(unsigned int h, unsigned int l)
  {
    if ( is_finite_by_mask_double(h) )
      return false;
    return
      (( h & CGAL_MANTISSA_DOUBLE_MASK ) != 0) || (( l & 0xffffffff ) != 0);
  }

public:
  static inline bool is_finite( const double& dble) {
    double d = dble;
    IEEE_754_double* p = reinterpret_cast<IEEE_754_double*>(&d);
    return is_finite_by_mask_double( p->c.H );
  }

  static inline bool is_valid( const double& dble) {
    double d = dble;
    IEEE_754_double* p = reinterpret_cast<IEEE_754_double*>(&d);
    return ! ( is_nan_by_mask_double( p->c.H, p->c.L ));
  }

#else
  static inline bool is_valid(double d) { return (d == d); }

  static inline bool is_finite(double d) {
    return (d == d) && (is_valid(d-d));
  }
#endif

  // GCC is faster with std::fabs().
#ifdef __GNUG__
  static inline double abs(const double& d) {
    return CGAL_CLIB_STD::fabs(d);
  }
#endif
};

#else // CGAL_NEW_NT_TRAITS

template<> struct Number_type_traits<double> {
  typedef Tag_false  Has_gcd;
  typedef Tag_true   Has_division;
  typedef Tag_true   Has_sqrt;
};

inline
const double &
to_double(const double & d)
{ return d; }

inline 
std::pair<double,double>
to_interval(double d)
{ return std::make_pair(d,d);}

inline
double
sqrt(double d)
{ return std::sqrt(d); }

#ifdef __sgi

inline
bool is_finite(double d)
{
    switch (fp_class_d(d)) {
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
bool is_valid(double d)
{
    switch (fp_class_d(d)) {
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

#elif defined CGAL_CFG_IEEE_754_BUG

//#define CGAL_EXPONENT_DOUBLE_MASK   0x7ff00000
//#define CGAL_MANTISSA_DOUBLE_MASK   0x000fffff

inline
bool
is_finite_by_mask_double(unsigned int h)
{
  unsigned int e = h & CGAL_EXPONENT_DOUBLE_MASK;
  return ( ( e ^ CGAL_EXPONENT_DOUBLE_MASK ) != 0 );
}

inline
bool
is_nan_by_mask_double(unsigned int h, unsigned int l)
{
  if ( is_finite_by_mask_double(h) )
      return false;
  return (( h & CGAL_MANTISSA_DOUBLE_MASK ) != 0) || (( l & 0xffffffff ) != 0);
}

inline
bool
is_finite( const double& dble)
{
  double d = dble;
  IEEE_754_double* p = reinterpret_cast<IEEE_754_double*>(&d);
  return is_finite_by_mask_double( p->c.H );
}

inline
bool
is_valid( const double& dble)
{
  double d = dble;
  IEEE_754_double* p = reinterpret_cast<IEEE_754_double*>(&d);
  return ! ( is_nan_by_mask_double( p->c.H, p->c.L ));
}

#else

inline
bool
is_valid(double d)
{ return (d == d); }

inline
bool
is_finite(double d)
{ return (d == d) && (is_valid(d-d)); }

#endif

inline
io_Operator
io_tag(double)
{ return io_Operator(); }

// GCC is faster with std::fabs().
#ifdef __GNUG__
inline
double
abs(const double& d)
{ return CGAL_CLIB_STD::fabs(d); }
#endif

#endif

CGAL_END_NAMESPACE

#endif // CGAL_DOUBLE_H
