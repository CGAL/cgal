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
// It will go away when the overlap_action() will be definable dynamically.
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

// Third template parameter is either
// No_Filter_Cache (the default) or Interval_nt_advanced.
struct No_Filter_Cache {};
typedef Interval_nt_advanced Filter_Cache;

template <class CT>
inline
void
compute_cache (const CT &value, No_Filter_Cache &cache)
{ }

template <class CT>
inline
void
compute_cache (const CT &value, Interval_nt_advanced &inter)
{ inter = convert_to<Interval_nt_advanced>(value); }

template <class CT>
inline
Interval_nt_advanced
give_interval (const CT &value, const No_Filter_Cache &cache)
{ return convert_to<Interval_nt_advanced>(value); }

template <class CT>
inline
Interval_nt_advanced
give_interval (const CT &value, const Interval_nt_advanced &inter)
{ return inter; }

template <class CT, class ET, class Cache_t = No_Filter_Cache >
struct Filtered_exact
{
  typedef Filtered_exact<CT,ET,Cache_t> Fil;
  typedef Interval_nt_advanced IA;

  Cache_t cache;
  CT value;
// #ifdef CGAL_FILTER_USE_CACHE
  void recompute_cache() { compute_cache (value, cache); }
  // { cache = convert_to<IA>(value); }
// #else
  // void recompute_cache() {}
// #endif

  Filtered_exact () {}
  template <class NT>
  Filtered_exact (const NT & nt)
      : value(nt)  { recompute_cache(); }
  // The following one is used for Quotient<>.
  template <class NT>
  Filtered_exact (const NT & num, const NT & den)
      : value(num, den)   { recompute_cache(); }
  Filtered_exact (const Filtered_exact<CT,ET,Cache_t> & fil)
      : value(fil.value)  { recompute_cache(); }
  // Filtered_exact (const double & d)	: value(d)  {}
  // Filtered_exact (const CT & ct)	: value(ct) {}

  // The two conversion functions are provided by the global scope.
  // Note that we could "cache" interval_value() [will be done].
  IA interval() const { return give_interval(value, cache); }
#if 0

#ifdef CGAL_FILTER_USE_CACHE
  IA interval() const { return cache; }
#else
#ifndef CGAL_CFG_NO_EXPLICIT_TEMPLATE_FUNCTION_ARGUMENT_SPECIFICATION
  IA interval() const { return convert_to<IA>(value); }
#else
  IA interval() const { return convert_from_to(IA(), value); }
#endif // CGAL_CFG_NO_EXPLICIT_TEMPLATE_FUNCTION_ARGUMENT_SPECIFICATION
#endif

#endif

#ifndef CGAL_CFG_NO_EXPLICIT_TEMPLATE_FUNCTION_ARGUMENT_SPECIFICATION
  ET exact()    const { return convert_to<ET>(value); }
#else
  ET exact()    const { return convert_from_to(ET(), value); }
#endif // CGAL_CFG_NO_EXPLICIT_TEMPLATE_FUNCTION_ARGUMENT_SPECIFICATION

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
template <class CT, class ET, class Cache_t>
inline
Filtered_exact<CT,ET,Cache_t>
sqrt (const Filtered_exact<CT,ET,Cache_t>& fil)
{ return sqrt(fil.value); }

template <class CT, class ET, class Cache_t>
inline
Filtered_exact<CT,ET,Cache_t>
square (const Filtered_exact<CT,ET,Cache_t>& fil)
{ return square(fil.value); }
#endif // CGAL_DENY_INEXACT_OPERATIONS_ON_FILTER

template <class CT, class ET, class Cache_t>
inline
bool
is_valid (const Filtered_exact<CT,ET,Cache_t>& fil)
{ return is_valid(fil.value); }

template <class CT, class ET, class Cache_t>
inline
bool
is_finite (const Filtered_exact<CT,ET,Cache_t>& fil)
{ return is_finite(fil.value); }

template <class CT, class ET, class Cache_t>
inline
double
to_double (const Filtered_exact<CT,ET,Cache_t>& fil)
{ return CGAL::to_double(fil.value); }

template <class CT, class ET, class Cache_t>
inline
Sign
sign (const Filtered_exact<CT,ET,Cache_t>& fil)
{ return CGAL::sign(fil.value); }

template <class CT, class ET, class Cache_t>
inline
Comparison_result
compare (const Filtered_exact<CT,ET,Cache_t>& fil,
	 const Filtered_exact<CT,ET,Cache_t>& fil2)
{ return CGAL::compare(fil.value, fil2.value); }

template <class CT, class ET, class Cache_t>
inline
Filtered_exact<CT,ET,Cache_t>
abs (const Filtered_exact<CT,ET,Cache_t>& fil)
{ return abs(fil.value); }

template <class CT, class ET, class Cache_t>
inline
Filtered_exact<CT,ET,Cache_t>
min (const Filtered_exact<CT,ET,Cache_t>& fil,
     const Filtered_exact<CT,ET,Cache_t>& fil2)
{ return min(fil.value, fil2.value); }

template <class CT, class ET, class Cache_t>
inline
Filtered_exact<CT,ET,Cache_t>
max (const Filtered_exact<CT,ET,Cache_t>& fil,
     const Filtered_exact<CT,ET,Cache_t>& fil2)
{ return max(fil.value, fil2.value); }

template <class CT, class ET, class Cache_t>
inline
io_Operator
io_tag (const Filtered_exact<CT,ET,Cache_t> &fil)
{ return io_tag(fil.value); }

template <class CT, class ET, class Cache_t>
inline
Number_tag
number_type_tag (const Filtered_exact<CT,ET,Cache_t> &fil)
{ return number_type_tag(fil.value); }

template <class CT, class ET, class Cache_t>
inline
std::ostream &
operator<< (std::ostream& os, const Filtered_exact<CT,ET,Cache_t>& d)
{ return os << d.value; }

template <class CT, class ET, class Cache_t>
inline
std::istream &
operator>> (std::istream &is, Filtered_exact<CT,ET,Cache_t>& d)
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
