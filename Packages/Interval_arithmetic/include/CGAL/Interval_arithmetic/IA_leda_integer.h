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
// file          : include/CGAL/Interval_arithmetic/IA_leda_integer.h
// revision      : $Revision$
// revision_date : $Date$
// package       : Interval Arithmetic
// author(s)     : Sylvain Pion <Sylvain.Pion@sophia.inria.fr>
//
// coordinator   : INRIA Sophia-Antipolis (<Mariette.Yvinec@sophia.inria.fr>)
//
// ============================================================================

#ifndef CGAL_IA_LEDA_INTEGER_H
#define CGAL_IA_LEDA_INTEGER_H

CGAL_BEGIN_NAMESPACE

// We choose the lazy approach, which is good enough: we take the double
// approximation, which is guaranted 1 bit error max, and return an interval
// around this value.  To have something more precise would require access to
// LEDA integer's internal representation, which is not possible without
// modifying LEDA.

inline
Interval_nt_advanced
convert_from_to (const Interval_nt_advanced&, const leda_integer & z)
{
#ifdef CGAL_IA_DEBUG
    CGAL_warning(FPU_get_cw() == FPU_cw_up);
#endif
    FPU_set_cw(FPU_cw_near);
    double approx = CGAL::to_double(z);
    FPU_set_cw(FPU_cw_up);
    Interval_nt_advanced result = approx + CGAL_IA_SMALLEST;
#ifdef CGAL_IA_DEBUG
    FPU_set_cw(FPU_cw_near);
    CGAL_warning(leda_integer(result.inf()) <= z &&
		 leda_integer(result.sup()) >= z);
    FPU_set_cw(FPU_cw_up);
#endif
    return result;
}

#ifndef CGAL_CFG_NO_EXPLICIT_TEMPLATE_FUNCTION_ARGUMENT_SPECIFICATION
template <>
struct converter<Interval_nt_advanced,leda_integer>
{
    static inline Interval_nt_advanced do_it (const leda_integer & z)
    {
	return convert_from_to(Interval_nt_advanced(), z);
    }
};
#endif // CGAL_CFG_NO_EXPLICIT_TEMPLATE_FUNCTION_ARGUMENT_SPECIFICATION


CGAL_END_NAMESPACE

#endif	 // CGAL_IA_LEDA_INTEGER_H
