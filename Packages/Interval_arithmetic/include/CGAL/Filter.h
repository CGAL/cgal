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
// file          : include/CGAL/Filter.h
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

#ifndef CGAL_FILTER_H
#define CGAL_FILTER_H

// CT = construction type (filtered)
// ET = exact type, used for exact predicate evaluation
// (Interval_nt_advanced) = used for filtering.
//
// 2 functions must be provided for the whole thing to work for a particular
// instantiation:
// - CGAL_to_interval_nt_advanced CGAL_to_interval_nt(const CT &);
//     which gives an interval surely containing the CT value.
// - ET CGAL_to_exact_type<ET>(const CT &)
//     which converts _exactly_ the CT value to ET.

template <class CT, class ET>
class CGAL_Filtering
{
public:
  CT value;

  CGAL_Filtering () {}
  CGAL_Filtering (int i) : value(i)  {}
  CGAL_Filtering (CT ct) : value(ct) {}

  typedef CGAL_Filtering<CT,ET> Fil;

  Fil operator-()               const { return Fil(-value); }

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

  bool operator< (const Fil& fil) const { return value <  fil.value; }
  bool operator> (const Fil& fil) const { return value >  fil.value; }
  bool operator<=(const Fil& fil) const { return value <= fil.value; }
  bool operator>=(const Fil& fil) const { return value >= fil.value; }
  bool operator==(const Fil& fil) const { return value == fil.value; }
  bool operator!=(const Fil& fil) const { return value != fil.value; }
};

template <class CT, class ET>
inline bool CGAL_is_valid    (const CGAL_Filtering<CT,ET>& fil)
{ return CGAL_is_valid(fil.value); }

template <class CT, class ET>
inline bool CGAL_is_finite   (const CGAL_Filtering<CT,ET>& fil)
{ return CGAL_is_finite(fil.value); }

template <class CT, class ET>
inline double CGAL_to_double (const CGAL_Filtering<CT,ET>& fil)
{ return CGAL_to_double(fil.value); }

#ifndef CGAL_DENY_INEXACT_OPERATIONS_ON_FILTER
template <class CT, class ET>
inline CGAL_Filtering<CT,ET> sqrt (const CGAL_Filtering<CT,ET>& fil)
{ return sqrt(fil.value); }
#endif // CGAL_DENY_INEXACT_OPERATIONS_ON_FILTER

// template <class ET>
  // template <class CT>
    // ET CGAL_to_exact_type(const CT &);

// All exact types should reasonnably have a built-in exact conversion
// from doubles ?  If not, it will fail, and you have to provide it.
//
// It's bad to provide such a default, because it can be an inexact cast:
// ex: CGAL_Gmpz accepts it, but it's false !!!

template <class ET>
inline ET CGAL_to_exact_type (const double & d)
{ return d; }


#ifdef CGAL_PREDICATES_ON_FTC2_H
#include <CGAL/Filter/predicates_on_ftC2.h>
#endif

#endif // CGAL_FILTER_H
