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
// file          : include/CGAL/Interval_arithmetic.h
// revision      : $Revision$
// revision_date : $Date$
// package       : Interval Arithmetic
// author(s)     : Sylvain Pion <Sylvain.Pion@sophia.inria.fr>
//
// coordinator   : INRIA Sophia-Antipolis (<Mariette.Yvinec@sophia.inria.fr>)
//
// ============================================================================

#ifndef CGAL_INTERVAL_ARITHMETIC_H
#define CGAL_INTERVAL_ARITHMETIC_H

// This file contains the description of the two classes:
// - Interval_nt_advanced  (do the FPU rounding mode changes yourself)
// - Interval_nt           ("plug-in" version, derived from the other one)
//
// The differences are:
// - The second one is slower.
// - The first one supposes the rounding mode is set -> +infinity before
// nearly all operations, and might set it -> +infinity when leaving, whereas
// the second leaves the rounding -> nearest.
//
// Note: When rounding is towards +infinity, to make an operation rounded
// towards -infinity, it's enough to take the opposite of some of the operand,
// and the opposite of the result (see operator+, operator*,...).

#include <CGAL/basic.h>
#include <CGAL/Interval_arithmetic/_FPU.h>	// FPU rounding mode functions.

#if defined _MSC_VER || defined __CYGWIN__
extern "C" { double CGAL_ms_sqrt(double); }
#endif

CGAL_BEGIN_NAMESPACE

#ifdef CGAL_IA_NO_INLINE
struct Interval_nt_advanced;
static Interval_nt_advanced operator*(const Interval_nt_advanced&,
	                              const Interval_nt_advanced&);
static Interval_nt_advanced operator/(const Interval_nt_advanced&,
	                              const Interval_nt_advanced&);
#endif

struct Interval_nt_advanced
{
  typedef Interval_nt_advanced IA;
  struct unsafe_comparison {};		// Exception class.
  static unsigned number_of_failures;	// Counts the number of failures.

  friend inline IA operator+     (const IA &, const IA &);
  friend inline IA operator-     (const IA &, const IA &);
#ifdef CGAL_IA_NO_INLINE
  friend        IA operator*     (const IA &, const IA &);
  friend        IA operator/     (const IA &, const IA &);
#else
  friend inline IA operator*     (const IA &, const IA &);
  friend inline IA operator/     (const IA &, const IA &);
#endif
  friend inline IA operator||    (const IA &, const IA &);
  friend inline IA operator&&    (const IA &, const IA &);
  friend inline bool operator<   (const IA &, const IA &);
  friend inline bool operator>   (const IA &, const IA &);
  friend inline bool operator<=  (const IA &, const IA &);
  friend inline bool operator>=  (const IA &, const IA &);
  friend inline bool operator==  (const IA &, const IA &);
  friend inline bool operator!=  (const IA &, const IA &);
  friend inline IA min           (const IA &, const IA &);
  friend inline IA max           (const IA &, const IA &);
  friend inline IA sqrt          (const IA &);
  friend inline IA square        (const IA &);
  friend inline IA abs           (const IA &);
  friend inline double to_double (const IA &);
  friend inline bool is_valid    (const IA &);
  friend inline bool is_finite   (const IA &);
  friend inline Sign sign        (const IA &);
  friend inline Comparison_result compare (const IA &, const IA &);

private:
  static void overlap_action() // throw (unsafe_comparison)
  {
      number_of_failures++;
      throw unsafe_comparison();
  }

public:

  // The constructors.
  Interval_nt_advanced() {}

  // To stop constant propagation, I need these CGAL_IA_STOP_CPROP().
  // Ideally, these barriers should be placed just before the critical
  // operations.

  Interval_nt_advanced(const double d)
  {
      _inf = _sup = CGAL_IA_STOP_CPROP(d);
  }

  Interval_nt_advanced(const double i, const double s)
  {
#ifndef CGAL_LAZY_EXACT_NT_H
      CGAL_assertion_msg(i<=s,
	      " CGAL bug or variable used before being initialized");
#endif
      _inf = CGAL_IA_STOP_CPROP(i);
      _sup = CGAL_IA_STOP_CPROP(s);
  }

#if 1
  // The copy constructors/assignment: useless.
  // The default ones are ok, but these appear to be faster...
  Interval_nt_advanced(const IA & d)
      : _inf(d._inf), _sup(d._sup) {}

  IA & operator=(const IA & d)
  { _inf = d._inf; _sup = d._sup; return *this; }
#endif

  IA  operator-() const { return IA (-_sup, -_inf); }

  IA & operator+= (const IA &);
  IA & operator-= (const IA &);
  IA & operator*= (const IA &);
  IA & operator/= (const IA &);

  bool is_same (const IA & d) const
  { return _inf == d._inf && _sup == d._sup; }

  bool is_point() const
  { return _sup == _inf; }

  bool overlap (const IA & d) const
  { return !(d._inf > _sup || d._sup < _inf); }

  double inf() const { return _inf; }
  double sup() const { return _sup; }

protected:
  double _inf, _sup;	// "_inf" stores the lower bound, "_sup" the upper.
};

// Two useful constant intervals.
// Smallest interval strictly containing zero.
static const Interval_nt_advanced Interval_Smallest
             (-CGAL_IA_MIN_DOUBLE, CGAL_IA_MIN_DOUBLE);
// [-inf;+inf]
static const Interval_nt_advanced Interval_Largest (-HUGE_VAL, HUGE_VAL);

// I'll remove those macros once it's tested.
#define CGAL_IA_SMALLEST CGAL::Interval_Smallest
#define CGAL_IA_LARGEST  CGAL::Interval_Largest

inline
Interval_nt_advanced
operator+ (const Interval_nt_advanced & e, const Interval_nt_advanced & d)
{
    CGAL_expensive_assertion(FPU_empiric_test() == FPU_cw_up);
    return Interval_nt_advanced (-CGAL_IA_FORCE_TO_DOUBLE((-e._inf) - d._inf),
	                          CGAL_IA_FORCE_TO_DOUBLE( e._sup + d._sup));
}

inline
Interval_nt_advanced
operator- (const Interval_nt_advanced & e, const Interval_nt_advanced & d)
{
    CGAL_expensive_assertion(FPU_empiric_test() == FPU_cw_up);
    return Interval_nt_advanced (-CGAL_IA_FORCE_TO_DOUBLE(d._sup - e._inf),
	                          CGAL_IA_FORCE_TO_DOUBLE(e._sup - d._inf));
}

#ifdef CGAL_IA_NO_INLINE
static
#else
inline
#endif
Interval_nt_advanced
operator* (const Interval_nt_advanced & e, const Interval_nt_advanced & d)
{
  CGAL_expensive_assertion(FPU_empiric_test() == FPU_cw_up);
  if (e._inf>=0.0)					// e>=0
  {
      // d>=0     [_inf*d._inf; _sup*d._sup]
      // d<=0     [_sup*d._inf; _inf*d._sup]
      // d~=0     [_sup*d._inf; _sup*d._sup]
      double a = e._inf, b = e._sup;
    if (d._inf < 0.0)
    {
	a=b;
	if (d._sup < 0.0)
	    b=e._inf;
    }

    return Interval_nt_advanced(-CGAL_IA_FORCE_TO_DOUBLE(a*(-d._inf)),
	                         CGAL_IA_FORCE_TO_DOUBLE(b*d._sup));
  }
  else if (e._sup<=0.0)				// e<=0
  {
      // d>=0     [_inf*d._sup; _sup*d._inf]
      // d<=0     [_sup*d._sup; _inf*d._inf]
      // d~=0     [_inf*d._sup; _inf*d._inf]
      double a = e._sup, b = e._inf;
    if (d._inf < 0.0)
    {
	a=b;
	if (d._sup < 0.0)
	    b=e._sup;
    }
    return Interval_nt_advanced(-CGAL_IA_FORCE_TO_DOUBLE(b*(-d._sup)),
	                         CGAL_IA_FORCE_TO_DOUBLE(a*d._inf));
  }
  else						// 0 \in [_inf;_sup]
  {
    if (d._inf>=0.0)				// d>=0
      return Interval_nt_advanced (-CGAL_IA_FORCE_TO_DOUBLE((-e._inf)*d._sup),
	                            CGAL_IA_FORCE_TO_DOUBLE(e._sup*d._sup));
    if (d._sup<=0.0)				// d<=0
      return Interval_nt_advanced (-CGAL_IA_FORCE_TO_DOUBLE(e._sup*(-d._inf)),
	                            CGAL_IA_FORCE_TO_DOUBLE(e._inf*d._inf));
        					// 0 \in d
    double tmp1 = CGAL_IA_FORCE_TO_DOUBLE((-e._inf)*d._sup);
    double tmp2 = CGAL_IA_FORCE_TO_DOUBLE(e._sup*(-d._inf));
    double tmp3 = CGAL_IA_FORCE_TO_DOUBLE(e._inf*d._inf);
    double tmp4 = CGAL_IA_FORCE_TO_DOUBLE(e._sup*d._sup);
    return Interval_nt_advanced(-std::max(tmp1,tmp2), std::max(tmp3,tmp4));
  };
}

#ifdef CGAL_IA_NO_INLINE
static
#else
inline
#endif
Interval_nt_advanced
operator/ (const Interval_nt_advanced & e, const Interval_nt_advanced & d)
{
  CGAL_expensive_assertion(FPU_empiric_test() == FPU_cw_up);
  if (d._inf>0.0)				// d>0
  {
      // e>=0	[_inf/d._sup; _sup/d._inf]
      // e<=0	[_inf/d._inf; _sup/d._sup]
      // e~=0	[_inf/d._inf; _sup/d._inf]
      double a = d._sup, b = d._inf;
    if (e._inf<0.0)
    {
	a=b;
	if (e._sup<0.0)
	    b=d._sup;
    };
    return Interval_nt_advanced(-CGAL_IA_FORCE_TO_DOUBLE((-e._inf)/a),
	                         CGAL_IA_FORCE_TO_DOUBLE(e._sup/b));
  }
  else if (d._sup<0.0)			// d<0
  {
      // e>=0	[_sup/d._sup; _inf/d._inf]
      // e<=0	[_sup/d._inf; _inf/d._sup]
      // e~=0	[_sup/d._sup; _inf/d._sup]
      double a = d._sup, b = d._inf;
    if (e._inf<0.0)
    {
	b=a;
	if (e._sup<0.0)
	    a=d._inf;
    };
    return Interval_nt_advanced(-CGAL_IA_FORCE_TO_DOUBLE((-e._sup)/a),
	                         CGAL_IA_FORCE_TO_DOUBLE(e._inf/b));
  }
  else					// d~0
    return CGAL_IA_LARGEST; // IA (-HUGE_VAL, HUGE_VAL);
	   // We could do slightly better -> [0;HUGE_VAL] when d._sup==0,
	   // but is this worth ?
}

inline
Interval_nt_advanced
sqrt (const Interval_nt_advanced & d)
{
  // sqrt([+a,+b]) => [sqrt(+a);sqrt(+b)]
  // sqrt([-a,+b]) => [0;sqrt(+b)] => assumes roundoff error.
  // sqrt([-a,-b]) => [0;sqrt(-b)] => assumes user bug (unspecified result).
  FPU_set_cw(FPU_cw_down);
#if defined _MSC_VER || defined __CYGWIN__
  // sqrt(double) on M$ is buggy.
  double i = (d._inf>0.0) ? CGAL_IA_FORCE_TO_DOUBLE(CGAL_ms_sqrt(d._inf)) : 0.0;
  FPU_set_cw(FPU_cw_up);
  return Interval_nt_advanced(i, CGAL_IA_FORCE_TO_DOUBLE(CGAL_ms_sqrt(d._sup)));
#else
  double i = (d._inf>0.0) ? CGAL_IA_FORCE_TO_DOUBLE(std::sqrt(d._inf)) : 0.0;
  FPU_set_cw(FPU_cw_up);
  return Interval_nt_advanced(i, CGAL_IA_FORCE_TO_DOUBLE(std::sqrt(d._sup)));
#endif
}

inline
Interval_nt_advanced
square (const Interval_nt_advanced & d)
{
  CGAL_expensive_assertion(FPU_empiric_test() == FPU_cw_up);
  if (d._inf>=0.0)
      return Interval_nt_advanced(-CGAL_IA_FORCE_TO_DOUBLE(d._inf*(-d._inf)),
	     			   CGAL_IA_FORCE_TO_DOUBLE(d._sup*d._sup));
  if (d._sup<=0.0)
      return Interval_nt_advanced(-CGAL_IA_FORCE_TO_DOUBLE(d._sup*(-d._sup)),
	     			   CGAL_IA_FORCE_TO_DOUBLE(d._inf*d._inf));
  return Interval_nt_advanced(0.0,
	  CGAL_IA_FORCE_TO_DOUBLE(square(std::max(-d._inf, d._sup))));
}

inline
bool
operator< (const Interval_nt_advanced & e, const Interval_nt_advanced & d)
{
    if (e._sup  < d._inf) return true;
    if (e._inf >= d._sup) return false;
    Interval_nt_advanced::overlap_action();
    return false;
}

inline
bool
operator<= (const Interval_nt_advanced & e, const Interval_nt_advanced & d)
{
    if (e._sup <= d._inf) return true;
    if (e._inf >  d._sup) return false;
    Interval_nt_advanced::overlap_action();
    return false;
}

inline
bool
operator== (const Interval_nt_advanced & e, const Interval_nt_advanced & d)
{
    if (d._inf >  e._sup || d._sup  < e._inf) return false;
    if (d._inf == e._sup && d._sup == e._inf) return true;
    Interval_nt_advanced::overlap_action();
    return false;
}

inline
bool
operator> (const Interval_nt_advanced & e, const Interval_nt_advanced & d)
{ return d < e; }

inline
bool
operator>= (const Interval_nt_advanced & e, const Interval_nt_advanced & d)
{ return d <= e; }

inline
bool
operator!= (const Interval_nt_advanced & e, const Interval_nt_advanced & d)
{ return !(d == e); }

inline
Interval_nt_advanced &
Interval_nt_advanced::operator+= (const Interval_nt_advanced & d)
{ return *this = *this + d; }

inline
Interval_nt_advanced &
Interval_nt_advanced::operator-= (const Interval_nt_advanced & d)
{ return *this = *this - d; }

inline
Interval_nt_advanced &
Interval_nt_advanced::operator*= (const Interval_nt_advanced & d)
{ return *this = *this * d; }

inline
Interval_nt_advanced &
Interval_nt_advanced::operator/= (const Interval_nt_advanced & d)
{ return *this = *this / d; }

inline
double
to_double (const Interval_nt_advanced & d)
{
    return (d._sup + d._inf) * 0.5;
}

inline
bool
is_valid (const Interval_nt_advanced & d)
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
is_finite (const Interval_nt_advanced & d)
{ return is_finite(d._inf) && is_finite(d._sup); }

inline
Sign
sign (const Interval_nt_advanced & d)
{
  if (d._inf > 0.0) return POSITIVE;
  if (d._sup < 0.0) return NEGATIVE;
  if (d._inf == d._sup) return ZERO;
  Interval_nt_advanced::overlap_action();
  return ZERO;
}

inline
Comparison_result
compare (const Interval_nt_advanced & d, const Interval_nt_advanced & e)
{
  if (d._inf > e._sup) return LARGER;
  if (e._inf > d._sup) return SMALLER;
  if (e._inf == d._sup && d._inf == e._sup) return EQUAL;
  Interval_nt_advanced::overlap_action();
  return EQUAL;
}

inline
Interval_nt_advanced
abs (const Interval_nt_advanced & d)
{
  if (d._inf >= 0.0) return d;
  if (d._sup <= 0.0) return -d;
  return Interval_nt_advanced(0.0, std::max(-d._inf, d._sup));
}

inline
Interval_nt_advanced
min (const Interval_nt_advanced & d, const Interval_nt_advanced & e)
{
  return Interval_nt_advanced(std::min(d._inf, e._inf),
			      std::min(d._sup, e._sup));
}

inline
Interval_nt_advanced
max (const Interval_nt_advanced & d, const Interval_nt_advanced & e)
{
  return Interval_nt_advanced(std::max(d._inf, e._inf),
			      std::max(d._sup, e._sup));
}

// The (join, union, ||) operator.
inline
Interval_nt_advanced
operator|| (const Interval_nt_advanced & d, const Interval_nt_advanced & e)
{
    return Interval_nt_advanced(std::min(e._inf, d._inf),
	                        std::max(e._sup, d._sup));
}

// The (meet, intersection, &&) operator.  Valid if intervals overlap.
inline
Interval_nt_advanced
operator&& (const Interval_nt_advanced & d, const Interval_nt_advanced & e)
{
    return Interval_nt_advanced(std::max(e._inf, d._inf),
	                        std::min(e._sup, d._sup));
}


// The non-advanced class.

struct Interval_nt : public Interval_nt_advanced
{
  typedef Interval_nt IA;

  // Constructors are identical.
  Interval_nt()
      {}
  Interval_nt(const double d)
      : Interval_nt_advanced(d) {}
  Interval_nt(const double a, const double b)
      : Interval_nt_advanced(a,b) {}

  // Private constructor for casts. (remade public)
  Interval_nt(const Interval_nt_advanced &d)
      : Interval_nt_advanced(d) {}

  IA operator-() const 
    { return IA(-_sup, -_inf); }

  IA & operator+=(const IA &);
  IA & operator-=(const IA &);
  IA & operator*=(const IA &);
  IA & operator/=(const IA &);
};


inline
Interval_nt
operator+ (const Interval_nt & e, const Interval_nt & d)
{
  FPU_CW_t backup = FPU_get_and_set_cw(FPU_cw_up);
  Interval_nt tmp ( Interval_nt_advanced(e) + Interval_nt_advanced(d) );
  FPU_set_cw(backup);
  return tmp;
}

inline
Interval_nt
operator- (const Interval_nt & e, const Interval_nt & d)
{
  FPU_CW_t backup = FPU_get_and_set_cw(FPU_cw_up);
  Interval_nt tmp ( Interval_nt_advanced(e) - Interval_nt_advanced(d) );
  FPU_set_cw(backup);
  return tmp;
}

inline
Interval_nt
operator* (const Interval_nt & e, const Interval_nt & d)
{
  FPU_CW_t backup = FPU_get_and_set_cw(FPU_cw_up);
  Interval_nt tmp ( Interval_nt_advanced(e) * Interval_nt_advanced(d) );
  FPU_set_cw(backup);
  return tmp;
}

inline
Interval_nt
operator/ (const Interval_nt & e, const Interval_nt & d)
{
  FPU_CW_t backup = FPU_get_and_set_cw(FPU_cw_up);
  Interval_nt tmp ( Interval_nt_advanced(e) / Interval_nt_advanced(d) );
  FPU_set_cw(backup);
  return tmp;
}

inline
Interval_nt
sqrt (const Interval_nt & d)
{
  FPU_CW_t backup = FPU_get_cw();
  Interval_nt tmp = sqrt( (Interval_nt_advanced) d);
  FPU_set_cw(backup);
  return tmp;
}

inline
Interval_nt
square (const Interval_nt & d)
{
  FPU_CW_t backup = FPU_get_and_set_cw(FPU_cw_up);
  Interval_nt tmp = square( (Interval_nt_advanced) d);
  FPU_set_cw(backup);
  return tmp;
}

inline
Interval_nt &
Interval_nt::operator+= (const Interval_nt & d)
{ return *this = *this + d; }

inline
Interval_nt &
Interval_nt::operator-= (const Interval_nt & d)
{ return *this = *this - d; }

inline
Interval_nt &
Interval_nt::operator*= (const Interval_nt & d)
{ return *this = *this * d; }

inline
Interval_nt &
Interval_nt::operator/= (const Interval_nt & d)
{ return *this = *this / d; }

inline
Interval_nt
abs (const Interval_nt & d)
{ return abs( (Interval_nt_advanced) d); }

inline
Interval_nt
min (const Interval_nt & d, const Interval_nt & e)
{ return min( (Interval_nt_advanced) d, (Interval_nt_advanced) e); }

inline
Interval_nt
max (const Interval_nt & d, const Interval_nt & e)
{ return max( (Interval_nt_advanced) d, (Interval_nt_advanced) e); }

inline
Interval_nt
operator|| (const Interval_nt & d, const Interval_nt & e)
{ return ((Interval_nt_advanced) d) || (Interval_nt_advanced) e; }

inline
Interval_nt
operator&& (const Interval_nt & d, const Interval_nt & e)
{ return ((Interval_nt_advanced) d) && (Interval_nt_advanced) e; }


std::ostream & operator<< (std::ostream &, const Interval_nt_advanced &);
std::istream & operator>> (std::istream &, Interval_nt_advanced &);

// The undocumented tags.

inline
io_Operator
io_tag (const Interval_nt_advanced &)
{ return io_Operator(); }

inline
Number_tag
number_type_tag (const Interval_nt_advanced &)
{ return Number_tag(); }

CGAL_END_NAMESPACE

// Finally we deal with the convert_to<Interval_nt_advanced>(NT).
//
// For the builtin types (well, all those that can be casted to double
// exactly), the template in misc.h is enough.

#ifdef CGAL_GMPZ_H
#include <CGAL/Interval_arithmetic/IA_Gmpz.h>
#endif

#ifdef CGAL_BIGFLOAT_H
#include <CGAL/Interval_arithmetic/IA_leda_bigfloat.h>
#endif

#ifdef CGAL_INTEGER_H
#include <CGAL/Interval_arithmetic/IA_leda_integer.h>
#endif

#ifdef CGAL_REAL_H
#include <CGAL/Interval_arithmetic/IA_leda_real.h>
#endif

#ifdef CGAL_RATIONAL_H
#include <CGAL/Interval_arithmetic/IA_leda_rational.h>
#endif

#ifdef CGAL_FIXED_PRECISION_NT_H
#include <CGAL/Interval_arithmetic/IA_Fixed.h>
#endif

#ifdef CGAL_QUOTIENT_H
#include <CGAL/Interval_arithmetic/IA_Quotient.h>
#endif

#endif // CGAL_INTERVAL_ARITHMETIC_H
