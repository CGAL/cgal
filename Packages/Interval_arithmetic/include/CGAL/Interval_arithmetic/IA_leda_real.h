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
// file          : include/CGAL/Interval_arithmetic/IA_leda_real.h
// revision      : $Revision$
// revision_date : $Date$
// package       : Interval Arithmetic
// author(s)     : Sylvain Pion <Sylvain.Pion@sophia.inria.fr>
//
// coordinator   : INRIA Sophia-Antipolis (<Mariette.Yvinec@sophia.inria.fr>)
//
// ============================================================================

#ifndef CGAL_IA_LEDA_REAL_H
#define CGAL_IA_LEDA_REAL_H

CGAL_BEGIN_NAMESPACE

inline
Interval_nt_advanced
convert_from_to (const Interval_nt_advanced&, const leda_real & z)
{
    CGAL_expensive_assertion(FPU_empiric_test() == CGAL_FE_UPWARD);
    FPU_set_cw(CGAL_FE_TONEAREST);
    double approx = CGAL::to_double(z);
    double rel_error = z.get_double_error();
    FPU_set_cw(CGAL_FE_UPWARD);
    Interval_nt_advanced result = approx
	* ( Interval_nt_advanced(-rel_error,rel_error) + 1 );
    CGAL_expensive_assertion_code(FPU_set_cw(CGAL_FE_TONEAREST);)
    CGAL_expensive_assertion( leda_real(result.inf()) <= z &&
		              leda_real(result.sup()) >= z );
    CGAL_expensive_assertion_code(FPU_set_cw(CGAL_FE_UPWARD);)
    return result;
}

#ifndef CGAL_CFG_NO_EXPLICIT_TEMPLATE_FUNCTION_ARGUMENT_SPECIFICATION
template <>
struct converter<Interval_nt_advanced,leda_real>
{
    static inline Interval_nt_advanced do_it (const leda_real & z)
    {
	return convert_from_to(Interval_nt_advanced(), z);
    }
};
#endif // CGAL_CFG_NO_EXPLICIT_TEMPLATE_FUNCTION_ARGUMENT_SPECIFICATION

CGAL_END_NAMESPACE

#endif	 // CGAL_IA_LEDA_REAL_H
