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
// file          : gmpxx.h
// package       : Number_types
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Sylvain Pion
//
// coordinator   : MPI, Saarbruecken
// ======================================================================
 
#ifndef CGAL_GMPXX_H
#define CGAL_GMPXX_H

#include <CGAL/basic.h>
#include <gmpxx.h>

// This file gathers the necessary adaptors so that the following
// C++ number types that come with GMP can be used by CGAL :
// - mpz_class
// - mpq_class
// - mpf_class

CGAL_BEGIN_NAMESPACE

template <>
struct Number_type_traits<mpz_class> {
  typedef Tag_false Has_gcd;
  typedef Tag_false Has_division;
  typedef Tag_true  Has_sqrt;
};

template <>
struct Number_type_traits<mpq_class> {
  typedef Tag_false Has_gcd;
  typedef Tag_true  Has_division;
  typedef Tag_false Has_sqrt;
};

template <>
struct Number_type_traits<mpf_class> {
  typedef Tag_false Has_gcd;
  typedef Tag_true  Has_division;
  typedef Tag_true  Has_sqrt;
};

inline
mpz_class
sqrt(const mpz_class &e)
{
    return ::sqrt(e);
}

inline
mpf_class
sqrt(const mpf_class &e)
{
    return ::sqrt(e);
}

inline
double
to_double(const mpz_class & e)
{ return e.get_d(); }

inline
double
to_double(const mpq_class & e)
{ return e.get_d(); }

inline
double
to_double(const mpf_class & e)
{ return e.get_d(); }

inline
bool
is_finite(const mpz_class &)
{ return true; }

inline
bool
is_valid(const mpz_class &)
{ return true; }

inline
bool
is_finite(const mpq_class &)
{ return true; }

inline
bool
is_valid(const mpq_class &)
{ return true; }

inline
bool
is_finite(const mpf_class &)
{ return true; }

inline
bool
is_valid(const mpf_class &)
{ return true; }

inline
io_Operator
io_tag(const mpz_class &)
{ return io_Operator(); }

inline
io_Operator
io_tag(const mpq_class &)
{ return io_Operator(); }

inline
io_Operator
io_tag(const mpf_class &)
{ return io_Operator(); }

inline
Sign
sign(const mpz_class& e)
{ return (Sign) ::sgn(e); }

inline
Sign
sign(const mpq_class& e)
{ return (Sign) ::sgn(e); }

inline
Sign
sign(const mpf_class& e)
{ return (Sign) ::sgn(e); }

inline
Comparison_result
compare(const mpz_class& e1, const mpz_class& e2)
{
  return (Comparison_result) ::cmp(e1, e2);
}

inline
Comparison_result
compare(const mpq_class& e1, const mpq_class& e2)
{
  return (Comparison_result) ::cmp(e1, e2);
}

inline
Comparison_result
compare(const mpf_class& e1, const mpf_class& e2)
{
  return (Comparison_result) ::cmp(e1, e2);
}

#if 0 // Unfinished
inline
Interval_base
to_interval (const mpq_class & e)
{
  Protect_FPU_rounding<true> P (CGAL_FE_TONEAREST);
  double approx = to_double(e);
  double rel_error = e.get_double_error();
  FPU_set_cw(CGAL_FE_UPWARD);
  return ( Interval_nt_advanced(-rel_error,rel_error) + 1 ) * approx;
}
#endif

namespace NTS {
  // These are necessary due to expression-templates.
  inline
  mpz_class
  abs(const mpz_class& x) { return ::abs(x); }

  inline
  mpq_class
  abs(const mpq_class& x) { return ::abs(x); }

  inline
  mpf_class
  abs(const mpf_class& x) { return ::abs(x); }
}

CGAL_END_NAMESPACE

#endif // CGAL_GMPXX_H
