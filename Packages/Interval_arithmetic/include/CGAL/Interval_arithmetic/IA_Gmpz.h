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
// coordinator   : INRIA Sophia-Antipolis (<Mariette.Yvinec@sophia.inria.fr>)
//
// ============================================================================

#ifndef CGAL_IA_GMPZ_H
#define CGAL_IA_GMPZ_H

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
	CGAL_expensive_assertion(FPU_empiric_test() == FPU_cw_up);
	FPU_set_cw(FPU_cw_near);
	double approx = CGAL::to_double(z);
	FPU_set_cw(FPU_cw_up);
	Interval_nt_advanced result = approx + CGAL_IA_SMALLEST;
	CGAL_expensive_assertion_code(FPU_set_cw(FPU_cw_near);)
	CGAL_expensive_assertion(Gmpz(result.inf()) <= z && Gmpz(result.sup()) >= z);
	CGAL_expensive_assertion_code(FPU_set_cw(FPU_cw_up);)
	return result;
}

#ifndef CGAL_CFG_NO_PARTIAL_TEMPLATE_FUNCTION_ARGUMENT_SPECIFICATION
template <>
struct converter<Interval_nt_advanced,Gmpz>
{
    static inline Interval_nt_advanced do_it (const Gmpz & z)
    {
	return convert_from_to(Interval_nt_advanced(), z);
    }
};
#endif // CGAL_CFG_NO_PARTIAL_TEMPLATE_FUNCTION_ARGUMENT_SPECIFICATION


CGAL_END_NAMESPACE

#endif	 // CGAL_IA_GMPZ_H
