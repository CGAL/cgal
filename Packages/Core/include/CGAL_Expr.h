/******************************************************************
 * Core Library Version 1.5, August 2002
 * Copyright (c) 1995-2002 Exact Computation Project
 * 
 * File: CGAL_Expr.h
 * Synopsis:
 *       Adds the necessary functionality so that CORE's Expr can 
 *       be used by CGAL.
 *
 * Written by
 *       Sylvain Pion <pion@cs.nyu.edu>
 *
 * WWW URL: http://cs.nyu.edu/exact/
 * Email: exact@cs.nyu.edu
 *
 * $Id$
 *****************************************************************/

// Adds the necessary functionality so that CORE::Expr can be used by CGAL.
// This file can also be found in CGAL itself.

#ifndef CORE_CGAL_EXPR_H
#define CORE_CGAL_EXPR_H

#include <CGAL/basic.h>
#define Level 4
#include <CORE.h>

CGAL_BEGIN_NAMESPACE

template <>
struct Number_type_traits<CORE::Expr> {
  typedef Tag_false Has_gcd;
  typedef Tag_true  Has_division;
  typedef Tag_true  Has_sqrt;
};

inline bool is_finite (const CORE::Expr&) { return true; }
inline bool is_valid (const CORE::Expr&) { return true; }
inline io_Operator io_tag (const CORE::Expr&) { return io_Operator(); }

inline double to_double (const CORE::Expr &e)
{ return e.doubleValue(); }

inline
CORE::Expr
sqrt(const CORE::Expr & e)
{ return CORE::sqrt(e); }

inline
Sign
sign(const CORE::Expr& e)
{ return (Sign) e.getSign(); }

inline
Comparison_result
compare(const CORE::Expr& e1, const CORE::Expr& e2)
{
  CORE::Expr c = e1-e2;
  return (c < 0) ? SMALLER : ((0 < c) ? LARGER : EQUAL);
}

#if 0 // unfinished
inline Interval_base to_interval (const CORE::Expr &e)
{
  bool this_need_to_be_filled;

  // Protect_FPU_rounding<true> P (CGAL_FE_TONEAREST);
  // Interval_nt_advanced approx(cl_double_approx(I));
  // FPU_set_cw(CGAL_FE_UPWARD);
  // return approx + Interval_base::Smallest;
}
#endif

// Specialized utilities :
namespace NTS {
  inline
  CORE::Expr
  sqrt( const CORE::Expr& e)
  { 
    return CORE::sqrt(e);
  }
} // namespace NTS

CGAL_END_NAMESPACE

#endif
