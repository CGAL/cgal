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
// file          : include/CGAL/Interval_arithmetic/IA_leda_bigfloat.h
// revision      : $Revision$
// revision_date : $Date$
// package       : Interval Arithmetic
// author(s)     : Sylvain Pion <Sylvain.Pion@sophia.inria.fr>
//
// coordinator   : INRIA Sophia-Antipolis (<Mariette.Yvinec@sophia.inria.fr>)
//
// ============================================================================

#ifndef CGAL_IA_LEDA_BIGFLOAT_H
#define CGAL_IA_LEDA_BIGFLOAT_H

// We choose the lazy approach, which is good enough: we take the double
// approximation, which is guaranted 1 bit error max(?), and return an
// interval around this value (+/- ulp).

inline
CGAL_Interval_nt_advanced
CGAL_convert_to (const leda_bigfloat &z, const CGAL_Interval_nt_advanced &)
{
#ifdef CGAL_IA_DEBUG
    CGAL_assertion(CGAL_FPU_get_rounding_mode() == CGAL_FPU_PLUS_INFINITY);
#endif
    CGAL_FPU_set_rounding_to_nearest();
    double approx = CGAL_to_double(z);
    CGAL_FPU_set_rounding_to_infinity();
    const CGAL_Interval_nt_advanced result = 
	CGAL_Interval_nt_advanced (approx) +
	CGAL_Interval_nt_advanced::smallest;
#ifdef CGAL_IA_DEBUG
    CGAL_FPU_set_rounding_to_nearest();
    CGAL_assertion(     leda_bigfloat(result.lower_bound()) <= z &&
			leda_bigfloat(result.upper_bound()) >= z);
    CGAL_FPU_set_rounding_to_infinity();
#endif
    return result;
}

#endif	 // CGAL_IA_LEDA_BIGFLOAT_H
