// Copyright (c) 1998-2004  Utrecht University (The Netherlands),
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
// Author(s)     : Sylvain Pion

#ifndef CGAL_FILTERED_EXACT_H
#define CGAL_FILTERED_EXACT_H
#define CGAL_ARITHMETIC_FILTER_H // for backward compatibility.

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
  : private Cache // To benefit from the empty base class optimization.
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

  void update_cache() { compute_cache (cache()); }

  const Cache & cache() const { return *this; }
  Cache & cache() { return *this; }

  // Private data member.
  CT    _value;

public:
#ifndef CGAL_NEW_NT_TRAITS
  typedef typename Number_type_traits<CT>::Has_gcd      Has_gcd;
  typedef typename Number_type_traits<CT>::Has_division Has_division;
  typedef typename Number_type_traits<CT>::Has_sqrt     Has_sqrt;
#endif

  Filtered_exact () {}
  Filtered_exact (const CT & ct)
      : _value(ct)  { update_cache(); }
  template <class NT>
  Filtered_exact (const NT & num, const NT & den) // For Quotient<>.
      : _value(num, den)   { update_cache(); }

  // The access functions.
  const CT & value() const { return _value; }
  IA interval() const { return give_interval(cache()); }
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


#ifdef CGAL_NEW_NT_TRAITS

template<class CT, class ET, bool Protected, class Cache>
struct Number_type_traits< Filtered_exact<CT,ET,Protected,Cache> >
  : public Number_type_traits<CT>
{
  typedef Filtered_exact<CT,ET,Protected,Cache> FENT;

#ifndef CGAL_DENY_INEXACT_OPERATIONS_ON_FILTER
  static inline FENT div (const FENT& fil1, const FENT& fil2)
  { 
    return CGAL::CGALi::checked_div(fil1.value(), fil2.value(), Has_gcd());
  }

  static inline FENT sqrt (const FENT& fil)
  { return CGAL::sqrt(fil.value()); }


  static inline FENT gcd (const FENT& fil1, const FENT& fil2)
  { 
    return CGAL::CGALi::checked_gcd(fil1.value(), fil2.value(), Has_gcd());
  }


  static inline FENT square (const FENT& fil)
  { return CGAL::square(fil.value()); }


#endif // CGAL_DENY_INEXACT_OPERATIONS_ON_FILTER

  static inline bool is_valid (const FENT& fil)
  { return CGAL::is_valid(fil.value()); }

  static inline bool is_finite (const FENT& fil)
  { return CGAL::is_finite(fil.value()); }

  static inline double to_double (const FENT& fil)
  { return CGAL::to_double(fil.value()); }

  static inline std::pair<double, double> to_interval (const FENT& fil)
  { return CGAL::to_interval(fil.value()); }

  static inline Sign sign (const FENT& fil)
  { return CGAL::sign(fil.value()); }


  static inline Comparison_result compare (const FENT& fil,
					   const FENT& fil2)
  { return CGAL::compare(fil.value(), fil2.value()); }


  static inline FENT abs (const FENT& fil)
  { return CGAL::abs(fil.value()); }


  static inline FENT min (const FENT& fil, const FENT& fil2)
  { return min(fil.value(), fil2.value()); }

  static inline FENT max (const FENT& fil, const FENT& fil2)
  { return max(fil.value(), fil2.value()); }

};

#else // CGAL_NEW_NT_TRAITS

// We forward the following functions to the CT value:
// sqrt, square, is_valid, is_finite, to_double, sign, compare, abs, min, max,
// div, gcd, io_tag, operator>>, operator<<.

#ifndef CGAL_DENY_INEXACT_OPERATIONS_ON_FILTER


template < class CT, class ET, bool Protected, class Cache >
inline
Filtered_exact<CT,ET,Dynamic,Protected,Cache>
div (const Filtered_exact<CT,ET,Dynamic, Protected,Cache>& fil1,
     const Filtered_exact<CT,ET, Dynamic, Protected,Cache>& fil2)
{ 
  typedef typename Filtered_exact<CT,
                                  ET,Dynamic,
                                  Protected,
                                  Cache>::Has_gcd Has_gcd;
  return CGAL::CGALi::checked_div(fil1.value(), fil2.value(), Has_gcd()); 
}


template < class CT, class ET, class Type, bool Protected, class Cache >
inline
Filtered_exact<CT, ET, Type, Protected, Cache>
sqrt (const Filtered_exact<CT, ET, Type, Protected, Cache>& fil)
{ return CGAL::sqrt(fil.value()); }


template < class CT, class ET, bool Protected, class Cache >
inline
Filtered_exact<CT,ET,Dynamic,Protected,Cache>
gcd (const Filtered_exact<CT,ET,Dynamic, Protected,Cache>& fil1,
     const Filtered_exact<CT,ET,Dynamic, Protected,Cache>& fil2)
{ 
  typedef typename Filtered_exact<CT,
                                  ET,Dynamic,
			          Protected,
                                  Cache>::Has_gcd Has_gcd;
  return CGAL::CGALi::checked_gcd(fil1.value(), fil2.value(), Has_gcd());
}


template < class CT, class ET, bool Protected, class Cache >
inline
Filtered_exact<CT,ET,Dynamic,Protected,Cache>
square (const Filtered_exact<CT, ET, Dynamic, Protected, Cache>& fil)
{ return CGAL_NTS square(fil.value()); }


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

template < class CT, class ET, class Type, bool Protected, class Cache >
inline
std::pair<double, double>
to_interval (const Filtered_exact<CT, ET, Type, Protected, Cache>& fil)
{ return CGAL::to_interval(fil.value()); }


template < class CT, class ET, bool Protected, class Cache >
inline
Sign
sign (const Filtered_exact<CT, ET, Dynamic, Protected, Cache>& fil)
{ return CGAL_NTS sign(fil.value()); }


template < class CT, class ET, bool Protected, class Cache >
inline
Comparison_result
compare (const Filtered_exact<CT, ET, Dynamic, Protected, Cache>& fil,
	 const Filtered_exact<CT, ET, Dynamic, Protected, Cache>& fil2)
{ return CGAL_NTS compare(fil.value(), fil2.value()); }


template < class CT, class ET, bool Protected, class Cache >
inline
Filtered_exact<CT,ET,Dynamic,Protected,Cache>
abs (const Filtered_exact<CT, ET, Dynamic, Protected, Cache>& fil)
{ return CGAL_NTS abs(fil.value()); }


template < class CT, class ET, bool Protected, class Cache >
inline
Filtered_exact<CT,ET,Dynamic,Protected,Cache>
min (const Filtered_exact<CT, ET, Dynamic, Protected, Cache>& fil,
     const Filtered_exact<CT, ET, Dynamic, Protected, Cache>& fil2)
{ return min(fil.value(), fil2.value()); }

template < class CT, class ET, bool Protected, class Cache >
inline
Filtered_exact<CT,ET,Dynamic,Protected,Cache>
max (const Filtered_exact<CT, ET, Dynamic, Protected, Cache>& fil,
     const Filtered_exact<CT, ET, Dynamic, Protected, Cache>& fil2)
{ return max(fil.value(), fil2.value()); }

template < class CT, class ET, class Type, bool Protected, class Cache >
inline
io_Operator
io_tag (const Filtered_exact<CT, ET, Type, Protected, Cache> &fil)
{ return io_tag(fil.value()); }

#endif // CGAL_NEW_NT_TRAITS

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

#endif // CGAL_FILTERED_EXACT_H
