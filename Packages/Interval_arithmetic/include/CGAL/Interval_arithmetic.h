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
#include <CGAL/FPU.h>
#include <CGAL/Interval_base.h>

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

  // The advantage of non-member operators is that (double * IA) just works...
  // But is it really useful and wishable in CGAL ?
  IA operator+ (const IA &d) const
  {
    Protect_FPU_rounding<Protected> P;
    return IA(-CGAL_IA_SUB(-inf_, d.inf_), CGAL_IA_ADD(sup_, d.sup_));
  }

  IA operator- (const IA &d) const
  {
    Protect_FPU_rounding<Protected> P;
    return IA(-CGAL_IA_SUB(d.sup_, inf_), CGAL_IA_SUB(sup_, d.inf_));
  }

  IA operator* (const IA &) const;
  IA operator/ (const IA &) const;

  IA operator-() const { return IA (-sup_, -inf_); }

  IA & operator+= (const IA &d) { return *this = *this + d; }
  IA & operator-= (const IA &d) { return *this = *this - d; }
  IA & operator*= (const IA &d) { return *this = *this * d; }
  IA & operator/= (const IA &d) { return *this = *this / d; }

  // The (join, union, ||) operator.
  IA operator|| (const IA & d) const
  {
    return IA(std::min(inf_, d.inf_), std::max(sup_, d.sup_));
  }

  // The (meet, intersection, &&) operator.  Valid if intervals overlap.
  IA operator&& (const IA & d) const
  {
    return IA(std::max(inf_, d.inf_), std::min(sup_, d.sup_));
  }
};

#ifdef LONG_LONG // missing CGAL_ ?
inline
Interval_base
to_interval (const long long & z)
{
  Protect_FPU_rounding<true> P(CGAL_FE_TONEAREST);
  Interval_nt_advanced approx (double(z));
  FPU_set_cw(CGAL_FE_UPWARD);
  return approx + Interval_base::Smallest;
}
#endif

template <bool Protected>
#ifndef CGAL_IA_NO_INLINE
inline
#endif
Interval_nt<Protected>
Interval_nt<Protected>::operator* (const Interval_nt<Protected> & d) const
{
  Protect_FPU_rounding<Protected> P;
  if (inf_>=0.0)					// e>=0
  {
    // d>=0     [inf_*d.inf_; sup_*d.sup_]
    // d<=0     [sup_*d.inf_; inf_*d.sup_]
    // d~=0     [sup_*d.inf_; sup_*d.sup_]
    double a = inf_, b = sup_;
    if (d.inf_ < 0.0)
    {
	a=b;
	if (d.sup_ < 0.0)
	    b=inf_;
    }
    return IA(-CGAL_IA_MUL(a, -d.inf_), CGAL_IA_MUL(b, d.sup_));
  }
  else if (sup_<=0.0)				// e<=0
  {
    // d>=0     [inf_*d.sup_; sup_*d.inf_]
    // d<=0     [sup_*d.sup_; inf_*d.inf_]
    // d~=0     [inf_*d.sup_; inf_*d.inf_]
    double a = sup_, b = inf_;
    if (d.inf_ < 0.0)
    {
	a=b;
	if (d.sup_ < 0.0)
	    b=sup_;
    }
    return IA(-CGAL_IA_MUL(b, -d.sup_), CGAL_IA_MUL(a, d.inf_));
  }
  else						// 0 \in [inf_;sup_]
  {
    if (d.inf_>=0.0)				// d>=0
      return IA(-CGAL_IA_MUL(-inf_, d.sup_), CGAL_IA_MUL(sup_, d.sup_));
    if (d.sup_<=0.0)				// d<=0
      return IA(-CGAL_IA_MUL(sup_, -d.inf_), CGAL_IA_MUL(inf_, d.inf_));
        					// 0 \in d
    double tmp1 = CGAL_IA_MUL(-inf_, d.sup_);
    double tmp2 = CGAL_IA_MUL(sup_, -d.inf_);
    double tmp3 = CGAL_IA_MUL(inf_, d.inf_);
    double tmp4 = CGAL_IA_MUL(sup_, d.sup_);
    return IA(-std::max(tmp1,tmp2), std::max(tmp3,tmp4));
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
  if (d.inf_>0.0)				// d>0
  {
    // e>=0	[inf_/d.sup_; sup_/d.inf_]
    // e<=0	[inf_/d.inf_; sup_/d.sup_]
    // e~=0	[inf_/d.inf_; sup_/d.inf_]
    double a = d.sup_, b = d.inf_;
    if (inf_<0.0)
    {
	a=b;
	if (sup_<0.0)
	    b=d.sup_;
    };
    return IA(-CGAL_IA_DIV(-inf_, a), CGAL_IA_DIV(sup_, b));
  }
  else if (d.sup_<0.0)			// d<0
  {
    // e>=0	[sup_/d.sup_; inf_/d.inf_]
    // e<=0	[sup_/d.inf_; inf_/d.sup_]
    // e~=0	[sup_/d.sup_; inf_/d.sup_]
    double a = d.sup_, b = d.inf_;
    if (inf_<0.0)
    {
	b=a;
	if (sup_<0.0)
	    a=d.inf_;
    };
    return IA(-CGAL_IA_DIV(-sup_, a), CGAL_IA_DIV(inf_, b));
  }
  else					// d~0
    return IA::Largest;
	   // We could do slightly better -> [0;HUGE_VAL] when d.sup_==0,
	   // but is this worth ?
}

#if 0 // Do this for the next release, same for is_one()
bool is_zero(const NT &n)
{
  if (0 >  n.sup_ || 0  < n.inf_) return false;
  if (0 == n.sup_ && 0 == n.inf_) return true;
  n.overlap_action();
}
#endif

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
  double i = (d.inf_ > 0.0) ? CGAL_IA_SQRT(d.inf_) : 0.0;
  FPU_set_cw(CGAL_FE_UPWARD);
  return Interval_nt<Protected>(i, CGAL_IA_SQRT(d.sup_));
}

#ifndef CGAL_CFG_MATCHING_BUG_2
template <bool Protected>
inline
Interval_nt<Protected>
min (const Interval_nt<Protected> & d, const Interval_nt<Protected> & e)
{
  return Interval_nt<Protected>(std::min(d.inf_, e.inf_),
		                std::min(d.sup_, e.sup_));
}

template <bool Protected>
inline
Interval_nt<Protected>
max (const Interval_nt<Protected> & d, const Interval_nt<Protected> & e)
{
  return Interval_nt<Protected>(std::max(d.inf_, e.inf_),
		                std::max(d.sup_, e.sup_));
}

#else

inline
Interval_nt<false>
min (const Interval_nt<false> & d, const Interval_nt<false> & e)
{
  return Interval_nt<false>(min(d.inf_, e.inf_),min(d.sup_, e.sup_));
}

inline
Interval_nt<false>
max (const Interval_nt<false> & d, const Interval_nt<false> & e)
{
  return Interval_nt<false>(max(d.inf_, e.inf_),max(d.sup_, e.sup_));
}

inline
Interval_nt<true>
min (const Interval_nt<true> & d, const Interval_nt<true> & e)
{
  return Interval_nt<true>(min(d.inf_, e.inf_),min(d.sup_, e.sup_));
}

inline
Interval_nt<true>
max (const Interval_nt<true> & d, const Interval_nt<true> & e)
{
  return Interval_nt<true>(max(d.inf_, e.inf_),max(d.sup_, e.sup_));
}

#endif // CGAL_CFG_MATCHING_BUG_2

namespace NTS {

#ifndef CGAL_CFG_MATCHING_BUG_2
template <bool Protected>
inline
Interval_nt<Protected>
square (const Interval_nt<Protected> & d)
{
  Protect_FPU_rounding<Protected> P;
  if (d.inf_>=0.0)
      return Interval_nt<Protected>(-CGAL_IA_MUL(d.inf_, -d.inf_),
	                             CGAL_IA_MUL(d.sup_, d.sup_));
  if (d.sup_<=0.0)
      return Interval_nt<Protected>(-CGAL_IA_MUL(d.sup_, -d.sup_),
	     	                     CGAL_IA_MUL(d.inf_, d.inf_));
  return Interval_nt<Protected>(0.0, CGAL_IA_SQUARE(std::max(-d.inf_,d.sup_)));
}

template <bool Protected>
inline
Interval_nt<Protected>
abs (const Interval_nt<Protected> & d)
{
  if (d.inf_ >= 0.0) return d;
  if (d.sup_ <= 0.0) return -d;
  return Interval_nt<Protected>(0.0, std::max(-d.inf_, d.sup_));
}

template <bool Protected>
inline
Sign
sign (const Interval_nt<Protected> & d)
{
  if (d.inf_ > 0.0) return POSITIVE;
  if (d.sup_ < 0.0) return NEGATIVE;
  if (d.inf_ == d.sup_) return ZERO;
  Interval_nt<Protected>::overlap_action();
  return ZERO;
}

template <bool Protected>
inline
Comparison_result
compare (const Interval_nt<Protected> & d, const Interval_nt<Protected> & e)
{
  if (d.inf_ > e.sup_) return LARGER;
  if (e.inf_ > d.sup_) return SMALLER;
  if (e.inf_ == d.sup_ && d.inf_ == e.sup_) return EQUAL;
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
  if (d.inf_>=0.0)
      return Interval_nt<true>(-CGAL_IA_MUL(d.inf_, -d.inf_),
	                        CGAL_IA_MUL(d.sup_, d.sup_));
  if (d.sup_<=0.0)
      return Interval_nt<true>(-CGAL_IA_MUL(d.sup_, -d.sup_),
	                        CGAL_IA_MUL(d.inf_, d.inf_));
  return Interval_nt<true>(0.0, CGAL_IA_SQUARE(std::max(-d.inf_, d.sup_)));
}

inline
Interval_nt<true>
abs (const Interval_nt<true> & d)
{
  if (d.inf_ >= 0.0) return d;
  if (d.sup_ <= 0.0) return -d;
  return Interval_nt<true>(0.0, std::max(-d.inf_, d.sup_));
}

inline
Sign
sign (const Interval_nt<true> & d)
{
  if (d.inf_ > 0.0) return POSITIVE;
  if (d.sup_ < 0.0) return NEGATIVE;
  if (d.inf_ == d.sup_) return ZERO;
  Interval_nt<true>::overlap_action();
  return ZERO;
}

inline
Comparison_result
compare (const Interval_nt<true> & d, const Interval_nt<true> & e)
{
  if (d.inf_ > e.sup_) return LARGER;
  if (e.inf_ > d.sup_) return SMALLER;
  if (e.inf_ == d.sup_ && d.inf_ == e.sup_) return EQUAL;
  Interval_nt<true>::overlap_action();
  return EQUAL;
}

// Then we overload for false.
inline
Interval_nt<false>
square (const Interval_nt<false> & d)
{
  Protect_FPU_rounding<false> P;
  if (d.inf_>=0.0)
      return Interval_nt<false>(-CGAL_IA_MUL(d.inf_, -d.inf_),
	                         CGAL_IA_MUL(d.sup_, d.sup_));
  if (d.sup_<=0.0)
      return Interval_nt<false>(-CGAL_IA_MUL(d.sup_, -d.sup_),
	                         CGAL_IA_MUL(d.inf_, d.inf_));
  return Interval_nt<false>(0.0, CGAL_IA_SQUARE(std::max(-d.inf_, d.sup_)));
}

inline
Interval_nt<false>
abs (const Interval_nt<false> & d)
{
  if (d.inf_ >= 0.0) return d;
  if (d.sup_ <= 0.0) return -d;
  return Interval_nt<false>(0.0, std::max(-d.inf_, d.sup_));
}

inline
Sign
sign (const Interval_nt<false> & d)
{
  if (d.inf_ > 0.0) return POSITIVE;
  if (d.sup_ < 0.0) return NEGATIVE;
  if (d.inf_ == d.sup_) return ZERO;
  Interval_nt<false>::overlap_action();
  return ZERO;
}

inline
Comparison_result
compare (const Interval_nt<false> & d, const Interval_nt<false> & e)
{
  if (d.inf_ > e.sup_) return LARGER;
  if (e.inf_ > d.sup_) return SMALLER;
  if (e.inf_ == d.sup_ && d.inf_ == e.sup_) return EQUAL;
  Interval_nt<false>::overlap_action();
  return EQUAL;
}
#endif // CGAL_CFG_MATCHING_BUG_2

} // namespace NTS

typedef Interval_nt<false> Interval_nt_advanced;  // for back-compatibility

CGAL_END_NAMESPACE

#endif // CGAL_INTERVAL_ARITHMETIC_H
