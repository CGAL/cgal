// ============================================================================
//
// Copyright (c) 2000 The CGAL Consortium
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
// file          : include/CGAL/Interval_base.h
// revision      : $Revision$
// revision_date : $Date$
// package       : Interval Arithmetic
// author(s)     : Sylvain Pion
// coordinator   : INRIA Sophia-Antipolis (<Mariette.Yvinec@sophia.inria.fr>)
//
// ============================================================================

#ifndef CGAL_INTERVAL_BASE_H
#define CGAL_INTERVAL_BASE_H

// This file contains the description of the Interval_base class, which is used
// by the number types Interval_nt<>.

#include <CGAL/basic.h>
#include <utility>				// Relational operators.

CGAL_BEGIN_NAMESPACE

struct Interval_base
{
  typedef Interval_base IA;
  struct unsafe_comparison {};		// Exception class.
  static unsigned number_of_failures;	// Number of filter failures.
  static const IA Smallest, Largest;	// Useful constant intervals.

  Interval_base () {}

  Interval_base (const double d)
    : _inf(d), _sup(d) {}

  Interval_base (const double i, const double s)
    : _inf(i), _sup(s)
  {
    CGAL_assertion_msg(i<=s,
	      " Variable used before being initialized (or CGAL bug)");
  }

  static void overlap_action() // throw (unsafe_comparison)
  {
    number_of_failures++;
    throw unsafe_comparison();
  }

  bool operator<  (const IA &d) const
  {
    if (_sup  < d._inf) return true;
    if (_inf >= d._sup) return false;
    overlap_action();
    return false;
  }

  bool operator<= (const IA &d) const
  {
    if (_sup <= d._inf) return true;
    if (_inf >  d._sup) return false;
    overlap_action();
    return false;
  }

  bool operator>= (const IA &d) const
  {
    return d <= *this;
  }

  bool operator== (const IA &d) const
  {
    if (d._inf >  _sup || d._sup  < _inf) return false;
    if (d._inf == _sup && d._sup == _inf) return true;
    overlap_action();
    return false;
  }

  bool is_point() const
  {
    return _sup == _inf;
  }

  bool is_same (const IA & d) const
  {
    return _inf == d._inf && _sup == d._sup;
  }

  bool do_overlap (const IA & d) const
  {
    return !(d._inf > _sup || d._sup < _inf);
  }

  double inf() const { return _inf; }
  double sup() const { return _sup; }

// protected:
  double _inf, _sup;	// "_inf" stores the lower bound, "_sup" the upper.
};

inline
double
to_double (const Interval_base & d)
{
  return (d._sup + d._inf) * 0.5;
}

inline
bool
is_valid (const Interval_base & d)
{
#if defined _MSC_VER || defined __sgi || defined __BORLANDC__
  return is_valid(d._inf) && is_valid(d._sup) && d._inf <= d._sup;
#else
  // The 2 first is_valid() are implicitely done by the 3rd test ;-)
  return d._inf <= d._sup;
#endif
}

inline
bool
is_finite (const Interval_base & d)
{
  return is_finite(d._inf) && is_finite(d._sup);
}

inline
io_Operator
io_tag (const Interval_base &)
{
  return io_Operator();
}

inline
Number_tag
number_type_tag (const Interval_base &)
{
  return Number_tag();
}

std::ostream & operator<< (std::ostream &, const Interval_base &);
std::istream & operator>> (std::istream &, Interval_base &);

CGAL_END_NAMESPACE

#endif // CGAL_INTERVAL_BASE_H
