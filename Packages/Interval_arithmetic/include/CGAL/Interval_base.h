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
#include <CGAL/double.h>                        // is_finite(double)

#ifdef __GNUG__
#  define CGAL_IA_NEW_FILTERS   // VC++ is not ready.
#endif

CGAL_BEGIN_NAMESPACE

struct Interval_base
{
  typedef Interval_base IA;
  struct unsafe_comparison {};		// Exception class.
  static unsigned number_of_failures;	// Number of filter failures.
  static const IA Smallest, Largest;	// Useful constant intervals.

  Interval_base () {}

  Interval_base (const double d)
    : inf_(d), sup_(d) {}

  Interval_base (const double i, const double s)
    : inf_(i), sup_(s)
  {
      // VC++ should use instead : (i<=s) || !is_valid(i) || !is_valid(s)
      // Or should I use is_valid() ? or is_valid_or_nan() ?
    CGAL_assertion_msg(!(i>s),
	      " Variable used before being initialized (or CGAL bug)");
  }

  static void overlap_action() // throw (unsafe_comparison)
  {
    number_of_failures++;
    throw unsafe_comparison();
  }

  bool operator<  (const IA &d) const
  {
    if (sup_  < d.inf_) return true;
    if (inf_ >= d.sup_) return false;
    overlap_action();
    return false;
  }

  bool operator>  (const IA &d) const
  {
    return d < *this;
  }

  bool operator<= (const IA &d) const
  {
    if (sup_ <= d.inf_) return true;
    if (inf_ >  d.sup_) return false;
    overlap_action();
    return false;
  }

  bool operator>= (const IA &d) const
  {
    return d <= *this;
  }

  bool operator== (const IA &d) const
  {
    if (d.inf_ >  sup_ || d.sup_  < inf_) return false;
    if (d.inf_ == sup_ && d.sup_ == inf_) return true;
    overlap_action();
    return false;
  }

  bool operator!= (const IA &d) const
  {
    return !(*this == d);
  }

  bool is_point() const
  {
    return sup_ == inf_;
  }

  bool is_same (const IA & d) const
  {
    return inf_ == d.inf_ && sup_ == d.sup_;
  }

  bool do_overlap (const IA & d) const
  {
    return !(d.inf_ > sup_ || d.sup_ < inf_);
  }

  double inf() const { return inf_; }
  double sup() const { return sup_; }

// protected:
  double inf_, sup_;	// "inf_" stores the lower bound, "sup_" the upper.
};

inline
Interval_base
to_interval (const Interval_base & d)
{
  return d;
}

// We put all the to_interval() of the builtin types here, because of #include
// circular dependencies otherwise.
inline
Interval_base
to_interval (const double & d)
{
  return Interval_base(d);
}

inline
Interval_base
to_interval (const float & f)
{
  return Interval_base(double(f));
}

inline
Interval_base
to_interval (const int & i)
{
  return Interval_base(double(i));
}

inline
Interval_base
to_interval (const short & s)
{
  return Interval_base(double(s));
}

inline
Interval_base
to_interval (const long & l)
{
  // actually we would like to compare number of mantissa bits,
  // this seems to be a sufficient approximation
  CGAL_assertion( sizeof(double) > sizeof(long) );
  // need something else for 64 bits longs.
  return Interval_base(double(l));
}

// to_interval(long long) is in Interval_arithmetic.h

inline
double
to_double (const Interval_base & d)
{
  return (d.sup_ + d.inf_) * 0.5;
  // This may overflow...
}

inline
bool
is_valid (const Interval_base & d)
{
#if defined _MSC_VER || defined __sgi || defined __BORLANDC__
  return CGAL::is_valid(d.inf_) && CGAL::is_valid(d.sup_) && d.inf_ <= d.sup_;
#else
  // The 2 first is_valid() are implicitely done by the 3rd test ;-)
  return d.inf_ <= d.sup_;
#endif
}

inline
bool
is_finite (const Interval_base & d)
{
  return CGAL::is_finite(d.inf_) && CGAL::is_finite(d.sup_);
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
