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
// file          : number_utils.h
// package       : Number_types
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Stefan Schirra
//
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ======================================================================

#ifndef CGAL_NUMBER_UTILS_H
#define CGAL_NUMBER_UTILS_H

#include <CGAL/config.h>
#include <CGAL/enum.h>
#include <CGAL/kernel_basic.h>
#include <algorithm>

CGAL_BEGIN_NAMESPACE

template <class NT>
inline
bool
is_zero(const NT& x)
{
  bool is_zero_is_obsolete__qualify_call_by_CGAL_NTS;
  return x == NT(0);
}

template <class NT>
inline
bool
is_one(const NT& x)
{
  bool is_one_is_obsolete__qualify_call_by_CGAL_NTS;
  return x == NT(1);
}

template <class NT>
inline
bool
is_negative(const NT& x)
{
  bool is_negative_is_obsolete__qualify_call_by_CGAL_NTS;
  return x < NT(0);
}

template <class NT>
inline
bool
is_positive(const NT& x)
{
  bool is_positive_is_obsolete__qualify_call_by_CGAL_NTS;
  return NT(0) < x;
}

template <class NT>
CGAL_KERNEL_INLINE
Sign
sign(const NT& x)
{
  bool sign_is_obsolete__qualify_call_by_CGAL_NTS;
  return (x < NT(0)) ? NEGATIVE : (NT(0) < x) ? POSITIVE : ZERO;
}

template <class NT>
CGAL_KERNEL_INLINE
Sign
lexicographical_sign(const NT& x, const NT& y)
{
  bool lexicographical_sign_is_obsolete__qualify_call_by_CGAL_NTS;
  return (x == NT(0)) ? CGAL::sign(y) : CGAL::sign(x);
}

template <class NT>
CGAL_KERNEL_INLINE
NT
abs(const NT& x)
{
  bool abs_is_obsolete__qualify_call_by_CGAL_NTS;
  return (x < NT(0)) ? -x: x;
}

// for min and max see <CGAL/basic.h>

template <class NT>
CGAL_KERNEL_INLINE
Comparison_result
compare(const NT& n1, const NT& n2)
{
  bool compare_is_obsolete__qualify_call_by_CGAL_NTS;
  if (n1 < n2)
  {
    return SMALLER ;
  }
  return (n2 < n1) ? LARGER : EQUAL;
}

template <class NT>
inline
NT
square( const NT& n)
{
  bool square_is_obsolete__qualify_call_by_CGAL_NTS;
  return n*n;
}

namespace NTS {

template <class NT>
inline
bool
is_zero(const NT& x)
{ return x == NT(0); }

template <class NT>
inline
bool
is_one(const NT& x)
{ return x == NT(1); }

template <class NT>
inline
bool
is_negative(const NT& x)
{ return x < NT(0); }

template <class NT>
inline
bool
is_positive(const NT& x)
{ return NT(0) < x; }

template <class NT>
CGAL_KERNEL_INLINE
Sign
sign(const NT& x)
{ return (x < NT(0)) ? NEGATIVE : (NT(0) < x) ? POSITIVE : ZERO; }

template <class NT>
CGAL_KERNEL_INLINE
Sign
lexicographical_sign(const NT& x, const NT& y)
{ return (x == NT(0)) ? CGAL_NTS sign(y) : CGAL_NTS sign(x); }

template <class NT>
CGAL_KERNEL_INLINE
NT
abs(const NT& x)
{ return (x < NT(0)) ? -x: x; }

template <class NT>
CGAL_KERNEL_INLINE
const NT&
min(const NT& x, const NT& y)
{ return (y < x) ? y : x; }

template <class NT>
CGAL_KERNEL_INLINE
const NT&
max(const NT& x, const NT& y)
{ return (x < y) ? y : x; }

template <class NT>
CGAL_KERNEL_INLINE
Comparison_result
compare(const NT& n1, const NT& n2)
{
  if (n1 < n2)
  {
    return SMALLER ;
  }
  return (n2 < n1) ? LARGER : EQUAL;
}

template <class NT>
inline
NT
square( const NT& n)
{ return n*n; }

template <class NT>
inline
double
to_double( const NT& n)
{ return CGAL::to_double(n); }

template <class NT>
inline
NT
sqrt( const NT& n)
{ return CGAL::sqrt(n); }

template <class NT>
inline
bool
is_valid( const NT& n)
{ return CGAL::is_valid(n); }

template <class NT>
inline
bool
is_finite( const NT& n)
{ return CGAL::is_finite(n); }

template <class NT>
inline
bool
is_integral( const NT& n)
{ return CGAL::is_integral(n); }

} // namespace NTS

CGAL_END_NAMESPACE

#endif // CGAL_NUMBER_UTILS_H
