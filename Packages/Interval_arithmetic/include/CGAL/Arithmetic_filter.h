// ============================================================================
//
// Copyright (c) 1998,1999 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       :
// release_date  :
//
// file          : include/CGAL/Arithmetic_filter.h
// revision      : $Revision$
// revision_date : $Date$
// package       : Interval Arithmetic
// author(s)     : Sylvain Pion <Sylvain.Pion@sophia.inria.fr>
//
// coordinator   : INRIA Sophia-Antipolis (<Mariette.Yvinec@sophia.inria.fr>)
//
// ============================================================================

// This file contains a wrapper type for number types, that helps for
// specializing template predicates, to use interval arithmetic as a filter.

#ifndef CGAL_ARITHMETIC_FILTER_H
#define CGAL_ARITHMETIC_FILTER_H

#include <CGAL/basic.h>

// #include <iostream>
#include <CGAL/Interval_arithmetic.h>
#if 0 // The following are already included from Interval_arithmetic.h.
#include <CGAL/enum.h>  // Because we overload {sign,compare,abs,min,max}
#include <CGAL/IO/io_tags.h>            // For io_tag().
#include <CGAL/number_type_tags.h>      // For number_type_tag()
#endif

// Check that the filtering stuff will work...
#ifdef CGAL_IA_NO_EXCEPTION
#  warning CGAL_IA_NO_EXCEPTION is defined !
#endif

CGAL_BEGIN_NAMESPACE

// CT = construction type
// ET = exact type (used for exact predicate evaluation)
// (Interval_nt_advanced) = used for filtering.
//
// 2 conversion functions must be provided:
// - convert_to <Interval_nt_advanced> (CT)
//     which gives an interval SURELY containing the CT value.
// - convert_to <ET> (CT)
//     which converts EXACTLY the CT value to ET.

template <class CT, class ET>
struct Filtered_exact
{
  typedef Filtered_exact<CT,ET> Fil;
  typedef Interval_nt_advanced IA;

  CT value;
#ifdef CGAL_FILTER_USE_CACHE
  IA cache;
  void recompute_cache() { cache = convert_to<IA>(value); }
#else
  void recompute_cache() {}
#endif

  Filtered_exact () {}
  template <class NT>
  Filtered_exact (const NT & nt)
      : value(nt)  { recompute_cache(); }
  // The following one is used for Quotient<>.
  template <class NT>
  Filtered_exact (const NT & num, const NT & den)
      : value(num, den)   { recompute_cache(); }
  Filtered_exact (const Filtered_exact<CT,ET> & fil)
      : value(fil.value)  { recompute_cache(); }
  // Filtered_exact (const double & d)	: value(d)  {}
  // Filtered_exact (const CT & ct)	: value(ct) {}

  // The two conversion functions are provided by the global scope.
  // Note that we could "cache" interval_value() [will be done].
#ifdef CGAL_FILTER_USE_CACHE
  IA interval() const { return cache; }
#else
  IA interval() const { return convert_to<IA>(value); }
#endif
  ET exact()    const { return convert_to<ET>(value); }

  // This one should not be needed, at least for now.
  // CT stored_value()   const { return value; }

  // Check Stroustrup if it's ok for assignment/ctors.
  Fil& operator= (const Fil& fil)
    { value = fil.value; recompute_cache(); return *this; }

  Fil  operator- ()               const { return Fil(-value); }
  bool operator< (const Fil& fil) const { return value <  fil.value; }
  bool operator> (const Fil& fil) const { return value >  fil.value; }
  bool operator<=(const Fil& fil) const { return value <= fil.value; }
  bool operator>=(const Fil& fil) const { return value >= fil.value; }
  bool operator==(const Fil& fil) const { return value == fil.value; }
  bool operator!=(const Fil& fil) const { return value != fil.value; }

#ifndef CGAL_DENY_INEXACT_OPERATIONS_ON_FILTER
  Fil operator+(const Fil& fil) const { return Fil(value + fil.value); }
  Fil operator-(const Fil& fil) const { return Fil(value - fil.value); }
  Fil operator*(const Fil& fil) const { return Fil(value * fil.value); }
  Fil operator/(const Fil& fil) const { return Fil(value / fil.value); }

  Fil& operator+=(const Fil& fil)
    { value += fil.value; recompute_cache(); return *this; }
  Fil& operator-=(const Fil& fil)
    { value -= fil.value; recompute_cache(); return *this; }
  Fil& operator*=(const Fil& fil)
    { value *= fil.value; recompute_cache(); return *this; }
  Fil& operator/=(const Fil& fil)
    { value /= fil.value; recompute_cache(); return *this; }
#endif // CGAL_DENY_INEXACT_OPERATIONS_ON_FILTER
};


// We forward the following functions to the CT value:
// sqrt, square, is_valid, is_finite, to_double, sign, compare, abs, min, max,
// io_tag, number_type_tag, operator>>, operator<<.

#ifndef CGAL_DENY_INEXACT_OPERATIONS_ON_FILTER
template <class CT, class ET>
inline
Filtered_exact<CT,ET>
sqrt (const Filtered_exact<CT,ET>& fil)
{ return sqrt(fil.value); }

template <class CT, class ET>
inline
Filtered_exact<CT,ET>
square (const Filtered_exact<CT,ET>& fil)
{ return square(fil.value); }
#endif // CGAL_DENY_INEXACT_OPERATIONS_ON_FILTER

template <class CT, class ET>
inline
bool
is_valid (const Filtered_exact<CT,ET>& fil)
{ return is_valid(fil.value); }

template <class CT, class ET>
inline
bool
is_finite (const Filtered_exact<CT,ET>& fil)
{ return is_finite(fil.value); }

template <class CT, class ET>
inline
double
to_double (const Filtered_exact<CT,ET>& fil)
{ return CGAL::to_double(fil.value); }

template <class CT, class ET>
inline
Sign
sign (const Filtered_exact<CT,ET>& fil)
{ return sign(fil.value); }

template <class CT, class ET>
inline
Comparison_result
compare (	const Filtered_exact<CT,ET>& fil,
		const Filtered_exact<CT,ET>& fil2)
{ return compare(fil.value, fil2.value); }

template <class CT, class ET>
inline
Filtered_exact<CT,ET>
abs (const Filtered_exact<CT,ET>& fil)
{ return abs(fil.value); }

template <class CT, class ET>
inline
Filtered_exact<CT,ET>
min (	const Filtered_exact<CT,ET>& fil,
		const Filtered_exact<CT,ET>& fil2)
{ return min(fil.value, fil2.value); }

template <class CT, class ET>
inline
Filtered_exact<CT,ET>
max (	const Filtered_exact<CT,ET>& fil,
		const Filtered_exact<CT,ET>& fil2)
{ return max(fil.value, fil2.value); }

template <class CT, class ET>
inline
io_Operator
io_tag (const Filtered_exact<CT,ET> &fil)
{ return io_tag(fil.value); }

template <class CT, class ET>
inline
Number_tag
number_type_tag (const Filtered_exact<CT,ET> &fil)
{ return number_type_tag(fil.value); }

template <class CT, class ET>
inline
std::ostream &
operator<< (std::ostream& os, const Filtered_exact<CT,ET>& d)
{ return os << d.value; }

template <class CT, class ET>
inline
std::istream &
operator>> (std::istream &is, Filtered_exact<CT,ET>& d)
{ return is >> d.value; }

CGAL_END_NAMESPACE

//  Now conditionnaly include the overloaded predicates.

#if defined( CGAL_PREDICATES_ON_FTC2_H ) && \
   !defined( CGAL_ARITHMETIC_FILTER_PREDICATES_ON_FTC2_H )
#include <CGAL/Arithmetic_filter/predicates_on_ftC2.h>
#endif

#if defined( CGAL_PREDICATES_ON_FTC3_H ) && \
   !defined( CGAL_ARITHMETIC_FILTER_PREDICATES_ON_FTC3_H )
#include <CGAL/Arithmetic_filter/predicates_on_ftC3.h>
#endif

#if defined( CGAL_PREDICATES_ON_RTH2_H ) && \
   !defined( CGAL_ARITHMETIC_FILTER_PREDICATES_ON_RTH2_H )
#include <CGAL/Arithmetic_filter/predicates_on_rtH2.h>
#endif

#if defined( CGAL_PREDICATES_ON_RTH3_H ) && \
   !defined( CGAL_ARITHMETIC_FILTER_PREDICATES_ON_RTH3_H )
#include <CGAL/Arithmetic_filter/predicates_on_rtH3.h>
#endif

#endif // CGAL_ARITHMETIC_FILTER_H
