// ======================================================================
//
// Copyright (c) 1999 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
// 
// release       : 
// release_date  : 
// 
// file          : double.h
// package       : Number_types
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Geert-Jan Giezeman
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ======================================================================

#ifndef CGAL_DOUBLE_H
#define CGAL_DOUBLE_H

#include <CGAL/basic.h>
#include <CGAL/tags.h>
#include <cmath>
#if defined(_MSC_VER) || defined(__BORLANDC__) || \
    defined(CGAL_MASK_FINITE_VALID)
#  include <CGAL/IEEE_754_unions.h>
#endif
#ifdef __sgi
#  include <fp_class.h>
#endif

CGAL_BEGIN_NAMESPACE

inline
double
to_double(double d)
{ return d; }

inline
double
sqrt(double d)
{ return std::sqrt(d); }

inline
bool
is_integral (const double d)
{ return ceil(d) == d; }

inline
Number_tag
number_type_tag(double)
{ return Number_tag(); }

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

#elif defined(_MSC_VER) || defined(__BORLANDC__) || \
      defined(CGAL_MASK_FINITE_VALID)

#define CGAL_EXPONENT_DOUBLE_MASK   0x7ff00000
#define CGAL_MANTISSA_DOUBLE_MASK   0x000fffff

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

namespace NTS {
#ifndef CGAL_NUMBER_UTILS_H
template <class NT> NT abs(const NT &x);
#endif // CGAL_NUMBER_UTILS_H

CGAL_TEMPLATE_NULL
inline
double
abs(const double& d)
{ return CGAL_CLIB_STD::fabs(d); }

} // namespace NTS

CGAL_END_NAMESPACE

#endif // CGAL_DOUBLE_H
