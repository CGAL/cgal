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
// file          : include/CGAL/Interval_arithmetic.h
// revision      : $Revision$
// revision_date : $Date$
// package       : Interval Arithmetic
// author(s)     : Sylvain Pion
// coordinator   : INRIA Sophia-Antipolis (<Mariette.Yvinec@sophia.inria.fr>)
//
// ============================================================================

#ifndef CGAL_INTERVAL_ARITHMETIC_H
#define CGAL_INTERVAL_ARITHMETIC_H

// This file contains the description of the following classes:
// - Interval_nt<false>  It's a number type that needs the FPU rounding mode
//                       to be set to +inf.  It is also typedef'd to
//                       Interval_nt_advanced for backward compatibility.
// - Interval_nt<true>   Same but it does the rounding mode itself so you
//                       don't have to worry about it.  But it's slower.
//
// Note: When rounding is towards +infinity, to make an operation rounded
// towards -infinity, it's enough to take the opposite of some of the operand,
// and the opposite of the result (see operator+, operator*,...).

#include <CGAL/basic.h>
#include <CGAL/Interval_arithmetic/_FPU.h>
#include <CGAL/Interval_base.h>

// sqrt(double) on M$ is buggy.
#if defined _MSC_VER || defined __CYGWIN__
extern "C" { double CGAL_ms_sqrt(double); }
#define CGAL_IA_SQRT(d) CGAL_ms_sqrt(d)
#else
#define CGAL_IA_SQRT(d) CGAL_CLIB_STD::sqrt(d)
#endif

CGAL_BEGIN_NAMESPACE

// bool or Tag_true/false, I don't know yet.
// Note that later, other behaviours might arise, such as "nothrow" => bool is
// probably not the best choice.  But it's easy to have !Protected. We'll see !
template <bool Protected = true>
struct Interval_nt : public Interval_base
{
  typedef Interval_nt<Protected> IA;

  Interval_nt() {}

  Interval_nt(const double d)
	  : Interval_base(d) {}

  Interval_nt(const double i, const double s)
	  : Interval_base(i,s) {}

  Interval_nt(const Interval_base & d)
	  : Interval_base(d) {}

#if 1
  // The copy constructors/assignment: useless.
  // The default ones are ok, but these appear to be faster with GCC 2.95.
  Interval_nt(const IA & d)
      : Interval_base(d._inf, d._sup) {}

  IA & operator=(const IA & d)
  { _inf = d._inf; _sup = d._sup; return *this; }
#endif

  // The advantage of non-member operators is that (double * IA) just works...
  // But is it really useful and wishable in CGAL ?
  IA operator+ (const IA &d) const
  {
    Protect_FPU_rounding<Protected> P;
    return Interval_nt<Protected>(-CGAL_IA_FORCE_TO_DOUBLE((-_inf) - d._inf),
	                           CGAL_IA_FORCE_TO_DOUBLE( _sup + d._sup));
  }

  IA operator- (const IA &d) const
  {
    Protect_FPU_rounding<Protected> P;
    return Interval_nt<Protected>(-CGAL_IA_FORCE_TO_DOUBLE(d._sup - _inf),
	                           CGAL_IA_FORCE_TO_DOUBLE(_sup - d._inf));
  }

  IA operator* (const IA &) const;
  IA operator/ (const IA &) const;

  IA operator-() const { return IA (-_sup, -_inf); }

  IA & operator+= (const IA &d) { return *this = *this + d; }
  IA & operator-= (const IA &d) { return *this = *this - d; }
  IA & operator*= (const IA &d) { return *this = *this * d; }
  IA & operator/= (const IA &d) { return *this = *this / d; }

  // The (join, union, ||) operator.
  IA operator|| (const IA & d) const
  {
    return Interval_nt<Protected>(std::min(_inf, d._inf),
                                  std::max(_sup, d._sup));
  }

  // The (meet, intersection, &&) operator.  Valid if intervals overlap.
  IA operator&& (const IA & d) const
  {
    return Interval_nt<Protected>(std::max(_inf, d._inf),
		                  std::min(_sup, d._sup));
  }
};

template <bool Protected>
#ifndef CGAL_IA_NO_INLINE
inline
#endif
Interval_nt<Protected>
Interval_nt<Protected>::operator* (const Interval_nt<Protected> & d) const
{
  Protect_FPU_rounding<Protected> P;
  if (_inf>=0.0)					// e>=0
  {
    // d>=0     [_inf*d._inf; _sup*d._sup]
    // d<=0     [_sup*d._inf; _inf*d._sup]
    // d~=0     [_sup*d._inf; _sup*d._sup]
    double a = _inf, b = _sup;
    if (d._inf < 0.0)
    {
	a=b;
	if (d._sup < 0.0)
	    b=_inf;
    }
    return Interval_nt<Protected>(-CGAL_IA_FORCE_TO_DOUBLE(a*(-d._inf)),
	                           CGAL_IA_FORCE_TO_DOUBLE(b*d._sup));
  }
  else if (_sup<=0.0)				// e<=0
  {
    // d>=0     [_inf*d._sup; _sup*d._inf]
    // d<=0     [_sup*d._sup; _inf*d._inf]
    // d~=0     [_inf*d._sup; _inf*d._inf]
    double a = _sup, b = _inf;
    if (d._inf < 0.0)
    {
	a=b;
	if (d._sup < 0.0)
	    b=_sup;
    }
    return Interval_nt<Protected>(-CGAL_IA_FORCE_TO_DOUBLE(b*(-d._sup)),
	                           CGAL_IA_FORCE_TO_DOUBLE(a*d._inf));
  }
  else						// 0 \in [_inf;_sup]
  {
    if (d._inf>=0.0)				// d>=0
      return Interval_nt<Protected>(-CGAL_IA_FORCE_TO_DOUBLE((-_inf)*d._sup),
	                             CGAL_IA_FORCE_TO_DOUBLE(_sup*d._sup));
    if (d._sup<=0.0)				// d<=0
      return Interval_nt<Protected>(-CGAL_IA_FORCE_TO_DOUBLE(_sup*(-d._inf)),
	                             CGAL_IA_FORCE_TO_DOUBLE(_inf*d._inf));
        					// 0 \in d
    double tmp1 = CGAL_IA_FORCE_TO_DOUBLE((-_inf)*d._sup);
    double tmp2 = CGAL_IA_FORCE_TO_DOUBLE(_sup*(-d._inf));
    double tmp3 = CGAL_IA_FORCE_TO_DOUBLE(_inf*d._inf);
    double tmp4 = CGAL_IA_FORCE_TO_DOUBLE(_sup*d._sup);
    return Interval_nt<Protected>(-std::max(tmp1,tmp2), std::max(tmp3,tmp4));
  };
}

template <bool Protected>
#ifndef CGAL_IA_NO_INLINE
inline
#endif
Interval_nt<Protected>
Interval_nt<Protected>::operator/ (const Interval_nt<Protected> & d) const
{
  Protect_FPU_rounding<Protected> P;
  if (d._inf>0.0)				// d>0
  {
    // e>=0	[_inf/d._sup; _sup/d._inf]
    // e<=0	[_inf/d._inf; _sup/d._sup]
    // e~=0	[_inf/d._inf; _sup/d._inf]
    double a = d._sup, b = d._inf;
    if (_inf<0.0)
    {
	a=b;
	if (_sup<0.0)
	    b=d._sup;
    };
    return Interval_nt<Protected>(-CGAL_IA_FORCE_TO_DOUBLE((-_inf)/a),
	                           CGAL_IA_FORCE_TO_DOUBLE(_sup/b));
  }
  else if (d._sup<0.0)			// d<0
  {
    // e>=0	[_sup/d._sup; _inf/d._inf]
    // e<=0	[_sup/d._inf; _inf/d._sup]
    // e~=0	[_sup/d._sup; _inf/d._sup]
    double a = d._sup, b = d._inf;
    if (_inf<0.0)
    {
	b=a;
	if (_sup<0.0)
	    a=d._inf;
    };
    return Interval_nt<Protected>(-CGAL_IA_FORCE_TO_DOUBLE((-_sup)/a),
	                           CGAL_IA_FORCE_TO_DOUBLE(_inf/b));
  }
  else					// d~0
    return Interval_nt<Protected>::Largest;
	   // We could do slightly better -> [0;HUGE_VAL] when d._sup==0,
	   // but is this worth ?
}

template <bool Protected>
inline
Interval_nt<Protected>
sqrt (const Interval_nt<Protected> & d)
{
  Protect_FPU_rounding<Protected> P;  // not optimal here.
  // sqrt([+a,+b]) => [sqrt(+a);sqrt(+b)]
  // sqrt([-a,+b]) => [0;sqrt(+b)] => assumes roundoff error.
  // sqrt([-a,-b]) => [0;sqrt(-b)] => assumes user bug (unspecified result).
  FPU_set_cw(CGAL_FE_DOWNWARD);
  double i = (d._inf > 0.0) ? CGAL_IA_FORCE_TO_DOUBLE(CGAL_IA_SQRT(d._inf))
	                    : 0.0;
  FPU_set_cw(CGAL_FE_UPWARD);
  return Interval_nt<Protected>
	  (i, CGAL_IA_FORCE_TO_DOUBLE(CGAL_IA_SQRT(d._sup)));
}

template <bool Protected>
inline
Interval_nt<Protected>
min (const Interval_nt<Protected> & d, const Interval_nt<Protected> & e)
{
  return Interval_nt<Protected>(std::min(d._inf, e._inf),
		                std::min(d._sup, e._sup));
}

template <bool Protected>
inline
Interval_nt<Protected>
max (const Interval_nt<Protected> & d, const Interval_nt<Protected> & e)
{
  return Interval_nt<Protected>(std::max(d._inf, e._inf),
		                std::max(d._sup, e._sup));
}

namespace NTS {

#ifndef CGAL_CFG_MATCHING_BUG_2
template <bool Protected>
inline
Interval_nt<Protected>
square (const Interval_nt<Protected> & d)
{
  Protect_FPU_rounding<Protected> P;
  if (d._inf>=0.0)
      return Interval_nt<Protected>(-CGAL_IA_FORCE_TO_DOUBLE(d._inf*(-d._inf)),
	     			     CGAL_IA_FORCE_TO_DOUBLE(d._sup*d._sup));
  if (d._sup<=0.0)
      return Interval_nt<Protected>(-CGAL_IA_FORCE_TO_DOUBLE(d._sup*(-d._sup)),
	     			     CGAL_IA_FORCE_TO_DOUBLE(d._inf*d._inf));
  return Interval_nt<Protected>(0.0,
	  CGAL_IA_FORCE_TO_DOUBLE(CGAL_NTS square(std::max(-d._inf, d._sup))));
}

template <bool Protected>
inline
Interval_nt<Protected>
abs (const Interval_nt<Protected> & d)
{
  if (d._inf >= 0.0) return d;
  if (d._sup <= 0.0) return -d;
  return Interval_nt<Protected>(0.0, std::max(-d._inf, d._sup));
}

template <bool Protected>
inline
Sign
sign (const Interval_nt<Protected> & d)
{
  if (d._inf > 0.0) return POSITIVE;
  if (d._sup < 0.0) return NEGATIVE;
  if (d._inf == d._sup) return ZERO;
  Interval_nt<Protected>::overlap_action();
  return ZERO;
}

template <bool Protected>
inline
Comparison_result
compare (const Interval_nt<Protected> & d, const Interval_nt<Protected> & e)
{
  if (d._inf > e._sup) return LARGER;
  if (e._inf > d._sup) return SMALLER;
  if (e._inf == d._sup && d._inf == e._sup) return EQUAL;
  Interval_nt<Protected>::overlap_action();
  return EQUAL;
}
#else // CGAL_CFG_MATCHING_BUG_2
// For crappy "compilers", we have to define complete overloaded functions.
// First we overload for true.
inline
Interval_nt<true>
square (const Interval_nt<true> & d)
{
  Protect_FPU_rounding<true> P;
  if (d._inf>=0.0)
      return Interval_nt<true>(-CGAL_IA_FORCE_TO_DOUBLE(d._inf*(-d._inf)),
	     			     CGAL_IA_FORCE_TO_DOUBLE(d._sup*d._sup));
  if (d._sup<=0.0)
      return Interval_nt<true>(-CGAL_IA_FORCE_TO_DOUBLE(d._sup*(-d._sup)),
	     			     CGAL_IA_FORCE_TO_DOUBLE(d._inf*d._inf));
  return Interval_nt<true>(0.0,
	  CGAL_IA_FORCE_TO_DOUBLE(CGAL_NTS square(std::max(-d._inf, d._sup))));
}

inline
Interval_nt<true>
abs (const Interval_nt<true> & d)
{
  if (d._inf >= 0.0) return d;
  if (d._sup <= 0.0) return -d;
  return Interval_nt<true>(0.0, std::max(-d._inf, d._sup));
}

inline
Sign
sign (const Interval_nt<true> & d)
{
  if (d._inf > 0.0) return POSITIVE;
  if (d._sup < 0.0) return NEGATIVE;
  if (d._inf == d._sup) return ZERO;
  Interval_nt<true>::overlap_action();
  return ZERO;
}

inline
Comparison_result
compare (const Interval_nt<true> & d, const Interval_nt<true> & e)
{
  if (d._inf > e._sup) return LARGER;
  if (e._inf > d._sup) return SMALLER;
  if (e._inf == d._sup && d._inf == e._sup) return EQUAL;
  Interval_nt<true>::overlap_action();
  return EQUAL;
}

// Then we overload for false.
inline
Interval_nt<false>
square (const Interval_nt<false> & d)
{
  Protect_FPU_rounding<false> P;
  if (d._inf>=0.0)
      return Interval_nt<false>(-CGAL_IA_FORCE_TO_DOUBLE(d._inf*(-d._inf)),
	     			     CGAL_IA_FORCE_TO_DOUBLE(d._sup*d._sup));
  if (d._sup<=0.0)
      return Interval_nt<false>(-CGAL_IA_FORCE_TO_DOUBLE(d._sup*(-d._sup)),
	     			     CGAL_IA_FORCE_TO_DOUBLE(d._inf*d._inf));
  return Interval_nt<false>(0.0,
	  CGAL_IA_FORCE_TO_DOUBLE(CGAL_NTS square(std::max(-d._inf, d._sup))));
}

inline
Interval_nt<false>
abs (const Interval_nt<false> & d)
{
  if (d._inf >= 0.0) return d;
  if (d._sup <= 0.0) return -d;
  return Interval_nt<false>(0.0, std::max(-d._inf, d._sup));
}

inline
Sign
sign (const Interval_nt<false> & d)
{
  if (d._inf > 0.0) return POSITIVE;
  if (d._sup < 0.0) return NEGATIVE;
  if (d._inf == d._sup) return ZERO;
  Interval_nt<false>::overlap_action();
  return ZERO;
}

inline
Comparison_result
compare (const Interval_nt<false> & d, const Interval_nt<false> & e)
{
  if (d._inf > e._sup) return LARGER;
  if (e._inf > d._sup) return SMALLER;
  if (e._inf == d._sup && d._inf == e._sup) return EQUAL;
  Interval_nt<false>::overlap_action();
  return EQUAL;
}
#endif // CGAL_CFG_MATCHING_BUG_2

} // namespace NTS

typedef Interval_nt<false> Interval_nt_advanced;  // for back-compatibility

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
