// ============================================================================
//
// Copyright (c) 1998 The CGAL Consortium
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

#include <iostream.h>	// Because we declare operator<< and >>.
#include <CGAL/enum.h>  // Because we overload CGAL_{sign,compare,abs,min,max}
#include <CGAL/IO/io_tags.h>            // For CGAL_io_Operator().
#include <CGAL/number_type_tags.h>      // For CGAL_number_type_tag()

// Check that the filtering stuff will work...
#ifdef CGAL_IA_NO_EXCEPTION
#  warning CGAL_IA_NO_EXCEPTION is defined !
#endif

// CT = construction type (filtered)
// ET = exact type, used for exact predicate evaluation
// (CGAL_Interval_nt_advanced) = used for filtering.
//
// 2 exact conversion functions must be provided:
// - CGAL_convert_to<CGAL_Interval_nt_advanced> (CT)
//     which gives an interval surely containing the CT value.
// - CGAL_convert_to<ET> (CT)
//     which converts _exactly_ the CT value to ET.

template <class CT, class ET>
struct CGAL_Filtered_exact
{
  CT value;

  CGAL_Filtered_exact () {}
  CGAL_Filtered_exact (const CGAL_Filtered_exact<CT,ET> & fil)
      : value(fil.value)  {}
  template <class NT>
  CGAL_Filtered_exact (const NT & nt) : value(nt)  {}
  // CGAL_Filtered_exact (const double & d)	: value(d)  {}
  // CGAL_Filtered_exact (const CT & ct)	: value(ct) {}

  typedef CGAL_Filtered_exact<CT,ET> Fil;

  // Check Stroustrup if it's ok for assignment/ctors.
  Fil& operator= (const Fil& fil) { value = fil.value; return *this; }

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

  Fil& operator+=(const Fil& fil) { value += fil.value; return *this; }
  Fil& operator-=(const Fil& fil) { value -= fil.value; return *this; }
  Fil& operator*=(const Fil& fil) { value *= fil.value; return *this; }
  Fil& operator/=(const Fil& fil) { value /= fil.value; return *this; }
#endif // CGAL_DENY_INEXACT_OPERATIONS_ON_FILTER
};


// We forward the following functions to the CT value:
// sqrt, CGAL_square, CGAL_is_valid, CGAL_is_finite, CGAL_to_double, CGAL_sign,
// CGAL_compare, CGAL_abs, CGAL_min, CGAL_max, CGAL_io_tag,
// CGAL_number_type_tag, operator>>, operator<<.

#ifndef CGAL_DENY_INEXACT_OPERATIONS_ON_FILTER
template <class CT, class ET>
inline
CGAL_Filtered_exact<CT,ET>
sqrt (const CGAL_Filtered_exact<CT,ET>& fil)
{ return sqrt(fil.value); }

template <class CT, class ET>
inline
CGAL_Filtered_exact<CT,ET>
CGAL_square (const CGAL_Filtered_exact<CT,ET>& fil)
{ return CGAL_square(fil.value); }
#endif // CGAL_DENY_INEXACT_OPERATIONS_ON_FILTER

template <class CT, class ET>
inline
bool
CGAL_is_valid (const CGAL_Filtered_exact<CT,ET>& fil)
{ return CGAL_is_valid(fil.value); }

template <class CT, class ET>
inline
bool
CGAL_is_finite (const CGAL_Filtered_exact<CT,ET>& fil)
{ return CGAL_is_finite(fil.value); }

template <class CT, class ET>
inline
double
CGAL_to_double (const CGAL_Filtered_exact<CT,ET>& fil)
{ return CGAL_to_double(fil.value); }

template <class CT, class ET>
inline
CGAL_Sign
CGAL_sign (const CGAL_Filtered_exact<CT,ET>& fil)
{ return CGAL_Sign(fil.value); }

template <class CT, class ET>
inline
CGAL_Comparison_result
CGAL_compare (	const CGAL_Filtered_exact<CT,ET>& fil,
		const CGAL_Filtered_exact<CT,ET>& fil2)
{ return CGAL_compare(fil.value, fil2.value); }

template <class CT, class ET>
inline
CGAL_Filtered_exact<CT,ET>
CGAL_abs (const CGAL_Filtered_exact<CT,ET>& fil)
{ return CGAL_abs(fil.value); }

template <class CT, class ET>
inline
CGAL_Filtered_exact<CT,ET>
CGAL_min (	const CGAL_Filtered_exact<CT,ET>& fil,
		const CGAL_Filtered_exact<CT,ET>& fil2)
{ return CGAL_min(fil.value, fil2.value); }

template <class CT, class ET>
inline
CGAL_Filtered_exact<CT,ET>
CGAL_max (	const CGAL_Filtered_exact<CT,ET>& fil,
		const CGAL_Filtered_exact<CT,ET>& fil2)
{ return CGAL_max(fil.value, fil2.value); }

template <class CT, class ET>
inline
CGAL_io_Operator
CGAL_io_tag (const CGAL_Filtered_exact<CT,ET> &fil)
{ return CGAL_io_tag(fil.value); }

template <class CT, class ET>
inline
CGAL_Number_tag
CGAL_number_type_tag (const CGAL_Filtered_exact<CT,ET> &fil)
{ return CGAL_number_type_tag(fil.value); }

template <class CT, class ET>
inline
ostream &
operator<< (ostream& os, const CGAL_Filtered_exact<CT,ET>& d)
{ return os << d.value; }

template <class CT, class ET>
inline
istream &
operator>> (istream &is, CGAL_Filtered_exact<CT,ET>& d)
{ return is >> d.value; }


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
