// ============================================================================
//
// Copyright (c) 1999 The CGAL Consortium
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
// file          : include/CGAL/Static_adaptatif_filter.h
// revision      : $Revision$
// revision_date : $Date$
// package       : Interval Arithmetic
// author(s)     : Sylvain Pion <Sylvain.Pion@sophia.inria.fr>
//
// coordinator   : INRIA Sophia-Antipolis (<Mariette.Yvinec@sophia.inria.fr>)
//
// ============================================================================

// This file contains a wrapper type for number types, that helps for
// specializing template predicates, to use the static adaptatif filter
// scheme.

#ifndef CGAL_STATIC_ADAPTATIF_FILTER_H
#define CGAL_STATIC_ADAPTATIF_FILTER_H

#include <CGAL/basic.h>
// #include <CGAL/Interval_arithmetic.h>

CGAL_BEGIN_NAMESPACE

// Currently this class is not template, but it might one day, just like
// Static_adaptatif_filter.
//
// Ideally, it should be integrated in Filtered_exact<> somehow
// (one more "enumerative" template parameter ?)

class Static_adaptatif_filter
{
  typedef Static_adaptatif_filter            Fil;
  // will get templated one day.
  // typedef typename Filtered_exact <double, leda_real> ET;
  typedef leda_real                          ET;
  typedef double			     CT;

  // Private data member.
  CT _value;

public:

  Static_adaptatif_filter () {}
  Static_adaptatif_filter (const Static_adaptatif_filter & fil)
      : _value(fil._value) {}
  Static_adaptatif_filter (const CT & d)
      : _value(d)  {}

  // The access functions.
  CT value() const { return _value; }
  ET exact() const { return ET(_value); }

  Fil  operator- ()               const { return Fil(-_value); }
  bool operator< (const Fil& fil) const { return _value <  fil._value; }
  bool operator> (const Fil& fil) const { return _value >  fil._value; }
  bool operator<=(const Fil& fil) const { return _value <= fil._value; }
  bool operator>=(const Fil& fil) const { return _value >= fil._value; }
  bool operator==(const Fil& fil) const { return _value == fil._value; }
  bool operator!=(const Fil& fil) const { return _value != fil._value; }

#ifndef CGAL_DENY_INEXACT_OPERATIONS_ON_FILTER
  Fil operator+(const Fil& fil) const { return Fil(_value + fil._value); }
  Fil operator-(const Fil& fil) const { return Fil(_value - fil._value); }
  Fil operator*(const Fil& fil) const { return Fil(_value * fil._value); }
  Fil operator/(const Fil& fil) const { return Fil(_value / fil._value); }

  Fil& operator+=(const Fil& fil)
    { _value += fil._value; return *this; }
  Fil& operator-=(const Fil& fil)
    { _value -= fil._value; return *this; }
  Fil& operator*=(const Fil& fil)
    { _value *= fil._value; return *this; }
  Fil& operator/=(const Fil& fil)
    { _value /= fil._value; return *this; }
#endif // CGAL_DENY_INEXACT_OPERATIONS_ON_FILTER
};


// We forward the following functions to the CT value:
// sqrt, square, is_valid, is_finite, to_double, sign, compare, abs, min, max,
// io_tag, number_type_tag, operator>>, operator<<.

#ifndef CGAL_DENY_INEXACT_OPERATIONS_ON_FILTER
inline
Static_adaptatif_filter
sqrt (const Static_adaptatif_filter& fil)
{ return std::sqrt(fil.value()); }

inline
Static_adaptatif_filter
square (const Static_adaptatif_filter& fil)
{ return square(fil.value()); }
#endif // CGAL_DENY_INEXACT_OPERATIONS_ON_FILTER

inline
bool
is_valid (const Static_adaptatif_filter& fil)
{ return is_valid(fil.value()); }

inline
bool
is_finite (const Static_adaptatif_filter& fil)
{ return is_finite(fil.value()); }

inline
double
to_double (const Static_adaptatif_filter& fil)
{ return CGAL::to_double(fil.value()); }

inline
Sign
sign (const Static_adaptatif_filter& fil)
{ return CGAL::sign(fil.value()); }

inline
Comparison_result
compare (const Static_adaptatif_filter& fil,
	 const Static_adaptatif_filter& fil2)
{ return CGAL::compare(fil.value(), fil2.value()); }

inline
Static_adaptatif_filter
abs (const Static_adaptatif_filter& fil)
{ return abs(fil.value()); }

inline
Static_adaptatif_filter
min (const Static_adaptatif_filter& fil,
     const Static_adaptatif_filter& fil2)
{ return std::min(fil.value(), fil2.value()); }

inline
Static_adaptatif_filter
max (const Static_adaptatif_filter& fil,
     const Static_adaptatif_filter& fil2)
{ return std::max(fil.value(), fil2.value()); }

inline
io_Operator
io_tag (const Static_adaptatif_filter &fil)
{ return io_tag(fil.value()); }

inline
Number_tag
number_type_tag (const Static_adaptatif_filter &fil)
{ return number_type_tag(fil.value()); }

inline
std::ostream &
operator<< (std::ostream& os, const Static_adaptatif_filter& d)
{ return os << d.value(); }

inline
std::istream &
operator>> (std::istream &is, Static_adaptatif_filter& d)
{
    typename Static_adaptatif_filter::CT e;
    is >> e;
    d = e; 
    return is;
}

CGAL_END_NAMESPACE

//  Now conditionnaly include the overloaded predicates.

#if defined( CGAL_PREDICATES_ON_FTC2_H ) && \
   !defined( CGAL_STATIC_ADAPTATIF_FILTER_PREDICATES_ON_FTC2_H )
#include <CGAL/Static_adaptatif_filter/predicates_on_ftC2.h>
#endif

#if defined( CGAL_PREDICATES_ON_FTC3_H ) && \
   !defined( CGAL_STATIC_ADAPTATIF_FILTER_PREDICATES_ON_FTC3_H )
#include <CGAL/Static_adaptatif_filter/predicates_on_ftC3.h>
#endif

#if defined( CGAL_PREDICATES_ON_RTH2_H ) && \
   !defined( CGAL_STATIC_ADAPTATIF_FILTER_PREDICATES_ON_RTH2_H )
#include <CGAL/Static_adaptatif_filter/predicates_on_rtH2.h>
#endif

#if defined( CGAL_PREDICATES_ON_RTH3_H ) && \
   !defined( CGAL_STATIC_ADAPTATIF_FILTER_PREDICATES_ON_RTH3_H )
#include <CGAL/Static_adaptatif_filter/predicates_on_rtH3.h>
#endif

#if defined( CGAL_REGULAR_TRIANGULATION_FTC2_H ) && \
   !defined( CGAL_STATIC_ADAPTATIF_FILTER_REGULAR_TRIANGULATION_FTC2_H )
#include <CGAL/Static_adaptatif_filter/predicates/Regular_triangulation_ftC2.h>
#endif

#if defined( CGAL_REGULAR_TRIANGULATION_RTH2_H ) && \
   !defined( CGAL_STATIC_ADAPTATIF_FILTER_REGULAR_TRIANGULATION_RTH2_H )
#include <CGAL/Static_adaptatif_filter/predicates/Regular_triangulation_rtH2.h>
#endif

#endif // CGAL_STATIC_ADAPTATIF_FILTER_H
