// ======================================================================
//
// Copyright (c) 2002 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
// 
// release       : 
// release_date  : 
// 
// file          : CORE_Expr.h
// package       : Number_types
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Sylvain Pion
//
// coordinator   : MPI, Saarbruecken
// ======================================================================
 
#ifndef CGAL_CORE_EXPR_H
#define CGAL_CORE_EXPR_H

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

inline
double
to_double(const CORE::Expr & e)
{ return e.doubleValue(); }

inline
CORE::Expr
sqrt(const CORE::Expr & e)
{ return CORE::sqrt(e); }

inline
bool
is_finite(const CORE::Expr &)
{ return true; }

inline
bool
is_valid(const CORE::Expr &)
{ return true; }

inline
io_Operator
io_tag(const CORE::Expr &)
{ return io_Operator(); }

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

#if 0 // Unfinished
inline
Interval_base
to_interval (const CORE::Expr & e)
{
  Protect_FPU_rounding<true> P (CGAL_FE_TONEAREST);
  double approx = to_double(e);
  double rel_error = e.get_double_error();
  FPU_set_cw(CGAL_FE_UPWARD);
  return ( Interval_nt_advanced(-rel_error,rel_error) + 1 ) * approx;
}
#endif

namespace NTS {
  inline
  CORE::Expr
  sqrt( const CORE::Expr& e)
  { 
    return CORE::sqrt(e);
  }
}

CGAL_END_NAMESPACE

#endif // CGAL_CORE_EXPR_H
