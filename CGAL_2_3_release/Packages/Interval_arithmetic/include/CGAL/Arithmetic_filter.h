// ============================================================================
//
// Copyright (c) 1998,1999,2000 The CGAL Consortium
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
// author(s)     : Sylvain Pion
// coordinator   : INRIA Sophia-Antipolis (<Mariette.Yvinec@sophia.inria.fr>)
//
// ============================================================================

#ifndef CGAL_ARITHMETIC_FILTER_H
#define CGAL_ARITHMETIC_FILTER_H

// This file contains a wrapper type for number types, that helps for
// specializing template predicates, to use interval arithmetic as a filter.

// Note: This stuff will one day be made obsolete by a filtered scheme at the
//       kernel level.

#include <CGAL/basic.h>
#include <CGAL/Interval_arithmetic.h>
#include <CGAL/Static_filter_error.h>
#include <CGAL/Restricted_double.h>
#include <CGAL/misc.h>

CGAL_BEGIN_NAMESPACE

// CT = construction type
// ET = exact type (used for exact predicate evaluation)
// Type = Static/Dynamic
// Protection = false/true (used to be Advanced/Protected)
// Cache = Filter_Cache/No_Filter_Cache
//
//   Static  + Protected : static adaptative
//   Static  + Advanced  : really static (needs some way to initialize bounds)
//   Dynamic + Protected : the current old one, based on intervals
//   Dynamic + Advanced  : rounding assumes +infty before entering predicates
//
// (Interval_nt_advanced) = used for filtering.
//
// 2 conversion functions must be provided:
// - to_interval(CT)
//     which gives an interval SURELY containing the CT value.
// - convert_to <ET> (CT)
//     which converts EXACTLY the CT value to ET.
//
// Let's add a .error() for when CT -> double is not exact.

// The user can initialize bounds using:
// Static_Filtered_orientationC2_6::new_bound(NEW_bound);

struct Static {};
struct Dynamic {};
struct No_Filter_Cache {};
typedef Interval_nt_advanced Filter_Cache;

template < class CT, class ET, class Type = Dynamic,
           bool Protection = true, class Cache = No_Filter_Cache >
class Filtered_exact
{
  typedef Filtered_exact<CT, ET, Type, Protection, Cache> Fil;
  typedef Interval_nt_advanced IA;

  // Cache managing functions.

  IA give_interval (const IA &inter) const
  {
      return inter;
  }
  IA give_interval (const No_Filter_Cache &) const
  {
      return Interval_nt_advanced(CGAL::to_interval(_value));
  }

  void compute_cache (const No_Filter_Cache &) const
  {}
  void compute_cache (IA &inter) const
  { inter = give_interval (No_Filter_Cache()); }

  void update_cache() { compute_cache (_cache); }

  // Private data members.

  CT _value;
  Cache _cache;

public:

  Filtered_exact () {}
  Filtered_exact (const Filtered_exact<CT, ET, Type, Protection, Cache> & fil)
      : _value(fil._value), _cache(fil._cache)  {}
  Filtered_exact (const CT & ct)
      : _value(ct)  { update_cache(); }
  template <class NT>
  Filtered_exact (const NT & num, const NT & den) // For Quotient<>.
      : _value(num, den)   { update_cache(); }

  // The access functions.
  CT value()	const { return _value; }
  IA interval() const { return give_interval(_cache); }
  ET exact()    const
  {
#ifndef CGAL_CFG_NO_PARTIAL_TEMPLATE_FUNCTION_ARGUMENT_SPECIFICATION
    return convert_to<ET>(_value);
#else
    return convert_from_to(ET(), _value);
#endif // CGAL_CFG_NO_PARTIAL_TEMPLATE_FUNCTION_ARGUMENT_SPECIFICATION
  }
  double to_double() const { return CGAL::to_double(_value); }
  Restricted_double dbl() const { return Restricted_double(to_double()); }

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
    { _value += fil._value; update_cache(); return *this; }
  Fil& operator-=(const Fil& fil)
    { _value -= fil._value; update_cache(); return *this; }
  Fil& operator*=(const Fil& fil)
    { _value *= fil._value; update_cache(); return *this; }
  Fil& operator/=(const Fil& fil)
    { _value /= fil._value; update_cache(); return *this; }
#endif // CGAL_DENY_INEXACT_OPERATIONS_ON_FILTER
};


// We forward the following functions to the CT value:
// sqrt, square, is_valid, is_finite, to_double, sign, compare, abs, min, max,
// io_tag, number_type_tag, operator>>, operator<<.

#ifndef CGAL_DENY_INEXACT_OPERATIONS_ON_FILTER
template <class CT, class ET, class Type, bool Protection, class Cache>
inline
Filtered_exact<CT,ET,Type,Protection,Cache>
sqrt (const Filtered_exact<CT,ET,Type,Protection,Cache>& fil)
{ return CGAL::sqrt(fil.value()); }

namespace NTS {

template <class CT, class ET, class Type, bool Protection, class Cache>
inline
Filtered_exact<CT,ET,Type,Protection,Cache>
square (const Filtered_exact<CT,ET,Type,Protection,Cache>& fil)
{ return CGAL_NTS square(fil.value()); }

} // namespace NTS
#endif // CGAL_DENY_INEXACT_OPERATIONS_ON_FILTER

template <class CT, class ET, class Type, bool Protection, class Cache>
inline
bool
is_valid (const Filtered_exact<CT,ET,Type,Protection,Cache>& fil)
{ return CGAL::is_valid(fil.value()); }

template <class CT, class ET, class Type, bool Protection, class Cache>
inline
bool
is_finite (const Filtered_exact<CT,ET,Type,Protection,Cache>& fil)
{ return CGAL::is_finite(fil.value()); }

template <class CT, class ET, class Type, bool Protection, class Cache>
inline
double
to_double (const Filtered_exact<CT,ET,Type,Protection,Cache>& fil)
{ return CGAL::to_double(fil.value()); }

namespace NTS {

template <class CT, class ET, class Type, bool Protection, class Cache>
inline
Sign
sign (const Filtered_exact<CT,ET,Type,Protection,Cache>& fil)
{ return CGAL_NTS sign(fil.value()); }

template <class CT, class ET, class Type, bool Protection, class Cache>
inline
Comparison_result
compare (const Filtered_exact<CT,ET,Type,Protection,Cache>& fil,
	 const Filtered_exact<CT,ET,Type,Protection,Cache>& fil2)
{ return CGAL_NTS compare(fil.value(), fil2.value()); }

template <class CT, class ET, class Type, bool Protection, class Cache>
inline
Filtered_exact<CT,ET,Type,Protection,Cache>
abs (const Filtered_exact<CT,ET,Type,Protection,Cache>& fil)
{ return CGAL_NTS abs(fil.value()); }

} // namespace NTS

template <class CT, class ET, class Type, bool Protection, class Cache>
inline
Filtered_exact<CT,ET,Type,Protection,Cache>
min (const Filtered_exact<CT,ET,Type,Protection,Cache>& fil,
     const Filtered_exact<CT,ET,Type,Protection,Cache>& fil2)
{ return min(fil.value(), fil2.value()); }

template <class CT, class ET, class Type, bool Protection, class Cache>
inline
Filtered_exact<CT,ET,Type,Protection,Cache>
max (const Filtered_exact<CT,ET,Type,Protection,Cache>& fil,
     const Filtered_exact<CT,ET,Type,Protection,Cache>& fil2)
{ return max(fil.value(), fil2.value()); }

template <class CT, class ET, class Type, bool Protection, class Cache>
inline
io_Operator
io_tag (const Filtered_exact<CT,ET,Type,Protection,Cache> &fil)
{ return io_tag(fil.value()); }

template <class CT, class ET, class Type, bool Protection, class Cache>
inline
Number_tag
number_type_tag (const Filtered_exact<CT,ET,Type,Protection,Cache> &fil)
{ return number_type_tag(fil.value()); }

template <class CT, class ET, class Type, bool Protection, class Cache>
inline
std::ostream &
operator<< (std::ostream& os,
	    const Filtered_exact<CT,ET,Type,Protection,Cache>& d)
{
    return os << d.value();
}

template <class CT, class ET, class Type, bool Protection, class Cache>
inline
std::istream &
operator>> (std::istream &is, Filtered_exact<CT,ET,Type,Protection,Cache>& d)
{
    CT e;
    is >> e;
    d = e; 
    return is;
}

CGAL_END_NAMESPACE

#include <CGAL/Arithmetic_filter/dispatch.h> // the overloaded predicates

#endif // CGAL_ARITHMETIC_FILTER_H
