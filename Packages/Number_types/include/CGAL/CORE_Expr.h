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
#include <CGAL/Number_type_traits.h>
#include <CGAL/Interval_arithmetic.h>

#include <utility>

#define CORE_LEVEL 4
#include <CORE/CORE.h>

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
{ return (Sign) e.sign(); }

inline
Comparison_result
compare(const CORE::Expr& e1, const CORE::Expr& e2)
{
  return Comparison_result(e1.cmp(e2));
}

// Should not be inline, but, well...
inline
std::pair<double,double>
to_interval (const CORE::Expr & e)
{
  std::pair<double,double> result;
  e.doubleInterval(result.first, result.second);
  return result;
}

CGAL_END_NAMESPACE

#endif // CGAL_CORE_EXPR_H
