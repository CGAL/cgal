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
// file          : include/CGAL/Interval_arithmetic/IA_Gmpz.h
// revision      : $Revision$
// revision_date : $Date$
// package       : Interval Arithmetic
// author(s)     : Sylvain Pion <Sylvain.Pion@sophia.inria.fr>
//
// ============================================================================

#ifndef CGAL_IA_GMPZ_H
#define CGAL_IA_GMPZ_H

#include <CGAL/Quotient.h> // Just for the converter double -> Quotient<Gmpz>.

CGAL_BEGIN_NAMESPACE

// We choose the lazy approach, which is good enough: we take the double
// approximation, which is guaranted 1 bit error max (when rounding to
// nearest), and return an interval around this value.
// It should be much faster to have a low level function especially designed
// for that using rounding to infinity.

inline
Interval_nt_advanced
convert_from_to (const Interval_nt_advanced&, const Gmpz & z)
{
	CGAL_expensive_assertion(FPU_empiric_test() == CGAL_FE_UPWARD);
	FPU_set_cw(CGAL_FE_TONEAREST);
	double approx = CGAL::to_double(z);
	FPU_set_cw(CGAL_FE_UPWARD);
	Interval_nt_advanced result = approx + Interval_nt_advanced::Smallest;
	CGAL_expensive_assertion_code(FPU_set_cw(CGAL_FE_TONEAREST);)
	CGAL_expensive_assertion(Gmpz(result.inf()) <= z &&
		                 Gmpz(result.sup()) >= z);
	CGAL_expensive_assertion_code(FPU_set_cw(CGAL_FE_UPWARD);)
	return result;
}

#ifndef CGAL_CFG_NO_EXPLICIT_TEMPLATE_FUNCTION_ARGUMENT_SPECIFICATION
template <>
struct converter<Interval_nt_advanced,Gmpz>
{
    static inline Interval_nt_advanced do_it (const Gmpz & z)
    {
	return convert_from_to(Interval_nt_advanced(), z);
    }
};
#endif // CGAL_CFG_NO_EXPLICIT_TEMPLATE_FUNCTION_ARGUMENT_SPECIFICATION


// Now we also have an exact converter from double to Quotient<Gmpz>, so that
// Filtered_exact<double, Quotient<Gmpz> > is useful.


// The following is an accessory function,
// which ideally should be moved to double.h.
// It tests if a double has an integral value.
// Result is unspecified for NaNs or Infs.
inline bool is_integral (const double d)
{
  return ceil(d) == d;
}

inline
Quotient<Gmpz>
convert_from_to (const Quotient<Gmpz>&, const double& d)
{ 
  // We multiply the value by 2, until it reaches an integral value.
  // Then it can be converted exactly to a Gmpz.
  // Note: it's not really optimized (it'll do 1000 iterations at worst).
  double num=d;
  double den=1.0;
  while ( ! CGAL::is_integral(num) )
  {
    den*=2.0;
    num*=2.0;
  }
  return Quotient<Gmpz>(Gmpz(num), Gmpz(den));
}

#ifndef CGAL_CFG_NO_EXPLICIT_TEMPLATE_FUNCTION_ARGUMENT_SPECIFICATION
template <>
struct converter<Quotient<Gmpz>,double>
{
  static inline Quotient<Gmpz> do_it (const double & z)
  { return convert_from_to(Quotient<Gmpz>(), z); }
};
#endif // CGAL_CFG_NO_EXPLICIT_TEMPLATE_FUNCTION_ARGUMENT_SPECIFICATION

CGAL_END_NAMESPACE

#endif // CGAL_IA_GMPZ_H
