// ============================================================================
//
// Copyright (c) 1998,1999,2000,2003 The CGAL Consortium
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
#include <CGAL/tags.h>
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
  typedef Interval_nt_advanced                            IA;

  // Cache managing functions.

  const IA & give_interval (const IA &inter) const
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
  {
    inter = give_interval (No_Filter_Cache());
  }

  void update_cache() { compute_cache (_cache); }

  // Private data members.

  CT    _value;
  Cache _cache;

public:

  typedef typename Number_type_traits<CT>::Has_gcd      Has_gcd;
  typedef typename Number_type_traits<CT>::Has_division Has_division;
  typedef typename Number_type_traits<CT>::Has_sqrt     Has_sqrt;

  Filtered_exact () {}
  Filtered_exact (const Filtered_exact<CT, ET, Type, Protection, Cache> & fil)
      : _value(fil._value), _cache(fil._cache)  {}
  Filtered_exact (const CT & ct)
      : _value(ct)  { update_cache(); }
  template <class NT>
  Filtered_exact (const NT & num, const NT & den) // For Quotient<>.
      : _value(num, den)   { update_cache(); }

  // The access functions.
  const CT & value() const { return _value; }
  IA interval() const { return give_interval(_cache); }
  ET exact()    const { return convert_to<ET>(_value); }

  double to_double() const { return CGAL::to_double(_value); }
  Restricted_double dbl() const { return Restricted_double(to_double()); }

  Fil  operator- ()  const { return Fil(-_value); }

#ifndef CGAL_DENY_INEXACT_OPERATIONS_ON_FILTER
  Fil& operator+=(const Fil& fil)
    { _value += fil._value; update_cache(); return *this; }
  Fil& operator-=(const Fil& fil)
    { _value -= fil._value; update_cache(); return *this; }
  Fil& operator*=(const Fil& fil)
    { _value *= fil._value; update_cache(); return *this; }
  Fil& operator/=(const Fil& fil)
    { _value /= fil._value; update_cache(); return *this; }
  Fil& operator%=(const Fil& fil)
    { _value %= fil._value; update_cache(); return *this; }
#endif // CGAL_DENY_INEXACT_OPERATIONS_ON_FILTER
};

#ifndef CGAL_DENY_INEXACT_OPERATIONS_ON_FILTER

template < class CT, class ET, class Type, bool Protected, class Cache >
inline
Filtered_exact<CT, ET, Type, Protected, Cache>
operator+(const Filtered_exact<CT, ET, Type, Protected, Cache>& a,
          const Filtered_exact<CT, ET, Type, Protected, Cache>& b)
{ return a.value() + b.value(); }

template < class CT, class ET, class Type, bool Protected, class Cache >
inline
Filtered_exact<CT, ET, Type, Protected, Cache>
operator-(const Filtered_exact<CT, ET, Type, Protected, Cache>& a,
          const Filtered_exact<CT, ET, Type, Protected, Cache>& b)
{ return a.value() - b.value(); }

template < class CT, class ET, class Type, bool Protected, class Cache >
inline
Filtered_exact<CT, ET, Type, Protected, Cache>
operator*(const Filtered_exact<CT, ET, Type, Protected, Cache>& a,
          const Filtered_exact<CT, ET, Type, Protected, Cache>& b)
{ return a.value() * b.value(); }

template < class CT, class ET, class Type, bool Protected, class Cache >
inline
Filtered_exact<CT, ET, Type, Protected, Cache>
operator/(const Filtered_exact<CT, ET, Type, Protected, Cache>& a,
          const Filtered_exact<CT, ET, Type, Protected, Cache>& b)
{ return a.value() / b.value(); }

// mixed operators
template < class CT, class ET, class Type, bool Protected, class Cache >
inline
Filtered_exact<CT, ET, Type, Protected, Cache>
operator+(const Filtered_exact<CT, ET, Type, Protected, Cache>& a,
          int b)
{ return a.value() + b; }

template < class CT, class ET, class Type, bool Protected, class Cache >
inline
Filtered_exact<CT, ET, Type, Protected, Cache>
operator-(const Filtered_exact<CT, ET, Type, Protected, Cache>& a,
          int b)
{ return a.value() - b; }

template < class CT, class ET, class Type, bool Protected, class Cache >
inline
Filtered_exact<CT, ET, Type, Protected, Cache>
operator*(const Filtered_exact<CT, ET, Type, Protected, Cache>& a,
          int b)
{ return a.value() * b; }

template < class CT, class ET, class Type, bool Protected, class Cache >
inline
Filtered_exact<CT, ET, Type, Protected, Cache>
operator/(const Filtered_exact<CT, ET, Type, Protected, Cache>& a,
          int b)
{ return a.value() / b; }

template < class CT, class ET, class Type, bool Protected, class Cache >
inline
Filtered_exact<CT, ET, Type, Protected, Cache>
operator+(int a,
          const Filtered_exact<CT, ET, Type, Protected, Cache>& b)
{ return a + b.value(); }

template < class CT, class ET, class Type, bool Protected, class Cache >
inline
Filtered_exact<CT, ET, Type, Protected, Cache>
operator-(int a,
          const Filtered_exact<CT, ET, Type, Protected, Cache>& b)
{ return a - b.value(); }

template < class CT, class ET, class Type, bool Protected, class Cache >
inline
Filtered_exact<CT, ET, Type, Protected, Cache>
operator*(int a,
          const Filtered_exact<CT, ET, Type, Protected, Cache>& b)
{ return a * b.value(); }

template < class CT, class ET, class Type, bool Protected, class Cache >
inline
Filtered_exact<CT, ET, Type, Protected, Cache>
operator/(int a,
          const Filtered_exact<CT, ET, Type, Protected, Cache>& b)
{ return a / b.value(); }

#endif // CGAL_DENY_INEXACT_OPERATIONS_ON_FILTER

template < class CT, class ET, class Type, bool Protected, class Cache >
inline
bool operator<(const Filtered_exact<CT, ET, Type, Protected, Cache>& a,
               const Filtered_exact<CT, ET, Type, Protected, Cache>& b)
{ return a.value() < b.value(); }

template < class CT, class ET, class Type, bool Protected, class Cache >
inline
bool operator>(const Filtered_exact<CT, ET, Type, Protected, Cache>& a,
               const Filtered_exact<CT, ET, Type, Protected, Cache>& b)
{ return a.value() > b.value(); }

template < class CT, class ET, class Type, bool Protected, class Cache >
inline
bool operator<=(const Filtered_exact<CT, ET, Type, Protected, Cache>& a,
                const Filtered_exact<CT, ET, Type, Protected, Cache>& b)
{ return a.value() <= b.value(); }

template < class CT, class ET, class Type, bool Protected, class Cache >
inline
bool operator>=(const Filtered_exact<CT, ET, Type, Protected, Cache>& a,
                const Filtered_exact<CT, ET, Type, Protected, Cache>& b)
{ return a.value() >= b.value(); }

template < class CT, class ET, class Type, bool Protected, class Cache >
inline
bool operator==(const Filtered_exact<CT, ET, Type, Protected, Cache>& a,
                const Filtered_exact<CT, ET, Type, Protected, Cache>& b)
{ return a.value() == b.value(); }

template < class CT, class ET, class Type, bool Protected, class Cache >
inline
bool operator!=(const Filtered_exact<CT, ET, Type, Protected, Cache>& a,
                const Filtered_exact<CT, ET, Type, Protected, Cache>& b)
{ return a.value() != b.value(); }


// mixed operators

template < class CT, class ET, class Type, bool Protected, class Cache >
inline
bool operator<(const Filtered_exact<CT, ET, Type, Protected, Cache>& a,
               int b)
{ return a.value() < b; }

template < class CT, class ET, class Type, bool Protected, class Cache >
inline
bool operator>(const Filtered_exact<CT, ET, Type, Protected, Cache>& a,
               int b)
{ return a.value() > b; }

template < class CT, class ET, class Type, bool Protected, class Cache >
inline
bool operator<=(const Filtered_exact<CT, ET, Type, Protected, Cache>& a,
                int b)
{ return a.value() <= b; }

template < class CT, class ET, class Type, bool Protected, class Cache >
inline
bool operator>=(const Filtered_exact<CT, ET, Type, Protected, Cache>& a,
                int b)
{ return a.value() >= b; }

template < class CT, class ET, class Type, bool Protected, class Cache >
inline
bool operator==(const Filtered_exact<CT, ET, Type, Protected, Cache>& a,
                int b)
{ return a.value() == b; }

template < class CT, class ET, class Type, bool Protected, class Cache >
inline
bool operator!=(const Filtered_exact<CT, ET, Type, Protected, Cache>& a,
                int b)
{ return a.value() != b; }

template < class CT, class ET, class Type, bool Protected, class Cache >
inline
bool operator<(int a,
               const Filtered_exact<CT, ET, Type, Protected, Cache>& b)
{ return a < b.value(); }

template < class CT, class ET, class Type, bool Protected, class Cache >
inline
bool operator>(int a,
               const Filtered_exact<CT, ET, Type, Protected, Cache>& b)
{ return a > b.value(); }

template < class CT, class ET, class Type, bool Protected, class Cache >
inline
bool operator<=(int a,
                const Filtered_exact<CT, ET, Type, Protected, Cache>& b)
{ return a <= b.value(); }

template < class CT, class ET, class Type, bool Protected, class Cache >
inline
bool operator>=(int a,
                const Filtered_exact<CT, ET, Type, Protected, Cache>& b)
{ return a >= b.value(); }

template < class CT, class ET, class Type, bool Protected, class Cache >
inline
bool operator==(int a,
                const Filtered_exact<CT, ET, Type, Protected, Cache>& b)
{ return a == b.value(); }

template < class CT, class ET, class Type, bool Protected, class Cache >
inline
bool operator!=(int a,
                const Filtered_exact<CT, ET, Type, Protected, Cache>& b)
{ return a != b.value(); }



namespace CGALi {

  // gcd and div make no sense for Filtered_exact with double.
  // As we have to provide a full specialization of gcd and div
  // for this number type we have to make it compile and make
  // it issue a warning at runtime.
  template <class NT>
  NT
  checked_gcd(const NT& n1, const NT& n2, Tag_true)
  {
    return CGAL_NTS gcd(n1,n2);
  }

  template <class NT>
  NT
  checked_gcd(const NT& n1, const NT& n2, Tag_false)
  {
    bool THIS_FUNCTION_SHOULD_ONLY_BE_INSTANTIATED_ON_MSVC_7_0;
    CGAL_assertion_msg(false, "numbertype has no gcd");
    return NT();
  }

  template <class NT>
  NT
  checked_div(const NT& n1, const NT& n2, Tag_true)
  {
    return CGAL_NTS div(n1,n2);
  }

  template <class NT>
  NT
  checked_div(const NT& n1, const NT& n2, Tag_false)
  {
    bool THIS_FUNCTION_SHOULD_ONLY_BE_INSTANTIATED_ON_MSVC_7_0;
    CGAL_assertion_msg(false, "numbertype has no div");
    return NT();
  }
} // namespace CGALi


// We forward the following functions to the CT value:
// sqrt, square, is_valid, is_finite, to_double, sign, compare, abs, min, max,
// div, gcd, io_tag, operator>>, operator<<.

#ifndef CGAL_DENY_INEXACT_OPERATIONS_ON_FILTER


#ifndef CGAL_CFG_MATCHING_BUG_2
template < class CGAL_IA_CT, class CGAL_IA_ET, bool CGAL_IA_PROTECTED,
           class CGAL_IA_CACHE >
inline
#else
static
#endif
Filtered_exact<CGAL_IA_CT,CGAL_IA_ET,Dynamic,CGAL_IA_PROTECTED,CGAL_IA_CACHE>
div (const Filtered_exact<CGAL_IA_CT,CGAL_IA_ET,Dynamic,
                          CGAL_IA_PROTECTED,CGAL_IA_CACHE>& fil1,
     const Filtered_exact<CGAL_IA_CT,CGAL_IA_ET, Dynamic,
                          CGAL_IA_PROTECTED,CGAL_IA_CACHE>& fil2)
{ return CGAL::CGALi::checked_div(fil1.value(), fil2.value(),
				  Filtered_exact<CGAL_IA_CT,CGAL_IA_ET,Dynamic,
				  CGAL_IA_PROTECTED,CGAL_IA_CACHE>::Has_gcd()); }


template < class CT, class ET, class Type, bool Protected, class Cache >
inline
Filtered_exact<CT, ET, Type, Protected, Cache>
sqrt (const Filtered_exact<CT, ET, Type, Protected, Cache>& fil)
{ return CGAL::sqrt(fil.value()); }



namespace NTS {

#ifndef CGAL_CFG_MATCHING_BUG_2
template < class CGAL_IA_CT, class CGAL_IA_ET, bool CGAL_IA_PROTECTED,
           class CGAL_IA_CACHE >
inline
#else
static
#endif
Filtered_exact<CGAL_IA_CT,CGAL_IA_ET,Dynamic,CGAL_IA_PROTECTED,CGAL_IA_CACHE>
gcd (const Filtered_exact<CGAL_IA_CT,CGAL_IA_ET,Dynamic,
                          CGAL_IA_PROTECTED,CGAL_IA_CACHE>& fil1,
     const Filtered_exact<CGAL_IA_CT,CGAL_IA_ET,Dynamic,
                          CGAL_IA_PROTECTED,CGAL_IA_CACHE>& fil2)
{ return CGAL::CGALi::checked_gcd(fil1.value(), fil2.value(), 
			    Filtered_exact<CGAL_IA_CT,CGAL_IA_ET,Dynamic,
				           CGAL_IA_PROTECTED,
				           CGAL_IA_CACHE>::Has_gcd()); }



#ifndef CGAL_CFG_MATCHING_BUG_2
template < class CGAL_IA_CT, class CGAL_IA_ET, bool CGAL_IA_PROTECTED,
           class CGAL_IA_CACHE >
inline
#else
static
#endif
Filtered_exact<CGAL_IA_CT,CGAL_IA_ET,Dynamic,CGAL_IA_PROTECTED,CGAL_IA_CACHE>
square (const Filtered_exact<CGAL_IA_CT, CGAL_IA_ET, Dynamic, 
                             CGAL_IA_PROTECTED, CGAL_IA_CACHE>& fil)
{ return CGAL_NTS square(fil.value()); }


} // namespace NTS
#endif // CGAL_DENY_INEXACT_OPERATIONS_ON_FILTER

template < class CT, class ET, class Type, bool Protected, class Cache >
inline
bool
is_valid (const Filtered_exact<CT, ET, Type, Protected, Cache>& fil)
{ return CGAL::is_valid(fil.value()); }

template < class CT, class ET, class Type, bool Protected, class Cache >
inline
bool
is_finite (const Filtered_exact<CT, ET, Type, Protected, Cache>& fil)
{ return CGAL::is_finite(fil.value()); }

template < class CT, class ET, class Type, bool Protected, class Cache >
inline
double
to_double (const Filtered_exact<CT, ET, Type, Protected, Cache>& fil)
{ return CGAL::to_double(fil.value()); }

namespace NTS {

#ifndef CGAL_CFG_MATCHING_BUG_2
template < class CGAL_IA_CT, class CGAL_IA_ET, bool CGAL_IA_PROTECTED,
           class CGAL_IA_CACHE >
inline
#else
static
#endif
Sign
sign (const Filtered_exact<CGAL_IA_CT, CGAL_IA_ET, Dynamic,
                           CGAL_IA_PROTECTED, CGAL_IA_CACHE>& fil)
{ return CGAL_NTS sign(fil.value()); }


#ifndef CGAL_CFG_MATCHING_BUG_2
template < class CGAL_IA_CT, class CGAL_IA_ET, bool CGAL_IA_PROTECTED,
           class CGAL_IA_CACHE >
inline
#else
static
#endif
Comparison_result
compare (const Filtered_exact<CGAL_IA_CT, CGAL_IA_ET, Dynamic,
                              CGAL_IA_PROTECTED, CGAL_IA_CACHE>& fil,
	 const Filtered_exact<CGAL_IA_CT, CGAL_IA_ET, Dynamic,
                              CGAL_IA_PROTECTED, CGAL_IA_CACHE>& fil2)
{ return CGAL_NTS compare(fil.value(), fil2.value()); }


#ifndef CGAL_CFG_MATCHING_BUG_2
template < class CGAL_IA_CT, class CGAL_IA_ET, bool CGAL_IA_PROTECTED,
           class CGAL_IA_CACHE >
inline
#else
static
#endif
Filtered_exact<CGAL_IA_CT,CGAL_IA_ET,Dynamic,CGAL_IA_PROTECTED,CGAL_IA_CACHE>
abs (const Filtered_exact<CGAL_IA_CT, CGAL_IA_ET, Dynamic,
                          CGAL_IA_PROTECTED, CGAL_IA_CACHE>& fil)
{ return CGAL_NTS abs(fil.value()); }

} // namespace NTS

#ifndef CGAL_CFG_MATCHING_BUG_2
template < class CGAL_IA_CT, class CGAL_IA_ET, bool CGAL_IA_PROTECTED,
           class CGAL_IA_CACHE >
inline
#else
static
#endif
Filtered_exact<CGAL_IA_CT,CGAL_IA_ET,Dynamic,CGAL_IA_PROTECTED,CGAL_IA_CACHE>
min (const Filtered_exact<CGAL_IA_CT, CGAL_IA_ET, Dynamic,
                          CGAL_IA_PROTECTED, CGAL_IA_CACHE>& fil,
     const Filtered_exact<CGAL_IA_CT, CGAL_IA_ET, Dynamic,
                          CGAL_IA_PROTECTED, CGAL_IA_CACHE>& fil2)
{ return min(fil.value(), fil2.value()); }

#ifndef CGAL_CFG_MATCHING_BUG_2
template < class CGAL_IA_CT, class CGAL_IA_ET, bool CGAL_IA_PROTECTED,
           class CGAL_IA_CACHE >
inline
#else
static
#endif
Filtered_exact<CGAL_IA_CT,CGAL_IA_ET,Dynamic,CGAL_IA_PROTECTED,CGAL_IA_CACHE>
max (const Filtered_exact<CGAL_IA_CT, CGAL_IA_ET, Dynamic,
                          CGAL_IA_PROTECTED, CGAL_IA_CACHE>& fil,
     const Filtered_exact<CGAL_IA_CT, CGAL_IA_ET, Dynamic,
                          CGAL_IA_PROTECTED, CGAL_IA_CACHE>& fil2)
{ return max(fil.value(), fil2.value()); }

template < class CT, class ET, class Type, bool Protected, class Cache >
inline
io_Operator
io_tag (const Filtered_exact<CT, ET, Type, Protected, Cache> &fil)
{ return io_tag(fil.value()); }

template < class CT, class ET, class Type, bool Protected, class Cache >
inline
std::ostream &
operator<< (std::ostream& os,
	    const Filtered_exact<CT, ET, Type, Protected, Cache>& d)
{ return os << d.value(); }

template < class CT, class ET, class Type, bool Protected, class Cache >
inline
std::istream &
operator>> (std::istream &is, 
	    Filtered_exact<CT, ET, Type, Protected, Cache>& d)
{
    CT e;
    is >> e;
    d = e; 
    return is;
}

CGAL_END_NAMESPACE

#include <CGAL/Arithmetic_filter/dispatch.h> // the overloaded predicates

#endif // CGAL_ARITHMETIC_FILTER_H
