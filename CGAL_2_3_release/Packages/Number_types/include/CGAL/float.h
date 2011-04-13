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
// file          : float.h
// package       : Number_types
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Geert-Jan Giezeman
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ======================================================================
 

#ifndef CGAL_FLOAT_H
#define CGAL_FLOAT_H

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
to_double(float f)
{ return static_cast<double>(f); }

inline
Number_tag
number_type_tag(float)
{ return Number_tag(); }

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
      defined(CGAL_MASK_FINITE_VALID)

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
