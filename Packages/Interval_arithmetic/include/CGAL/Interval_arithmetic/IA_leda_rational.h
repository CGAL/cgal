// ============================================================================
//
// Copyright (c) 1998 The CGAL Consortium
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
// file          : include/CGAL/Interval_arithmetic/IA_leda_rational.h
// revision      : $Revision$
// revision_date : $Date$
// package       : Interval Arithmetic
// author(s)     : Sylvain Pion <Sylvain.Pion@sophia.inria.fr>
//
// coordinator   : INRIA Sophia-Antipolis (<Mariette.Yvinec@sophia.inria.fr>)
//
// ============================================================================

#ifndef CGAL_IA_LEDA_RATIONAL_H
#define CGAL_IA_LEDA_RATIONAL_H

// For this one, I hope that adding 3 ulps will be enough for an exact
// conversion.  Since LEDA types (except real) don't give information on the
// precision of to_double(), we can't do much...

inline CGAL_Interval_nt_advanced CGAL_convert_to (const leda_rational &z)
{
#ifndef CGAL_NO_PRECONDITIONS
    CGAL_assertion(CGAL_FPU_get_rounding_mode() == CGAL_FPU_PLUS_INFINITY);
#endif
    CGAL_FPU_set_rounding_to_nearest();
    double approx = CGAL_to_double(z);
    CGAL_FPU_set_rounding_to_infinity();

    const CGAL_Interval_nt_advanced result =
	((CGAL_Interval_nt_advanced (approx)
	  + CGAL_Interval_nt_advanced::smallest)
	  + CGAL_Interval_nt_advanced::smallest)
	  + CGAL_Interval_nt_advanced::smallest;
// The following is bad because overflow is highly probable with rationals.
    // return CGAL_convert_to<CGAL_Interval_nt_advanced>(z.numerator())
	// /  CGAL_convert_to<CGAL_Interval_nt_advanced>(z.denominator());
#ifndef CGAL_NO_POSTCONDITIONS
    CGAL_FPU_set_rounding_to_nearest();
    CGAL_assertion( leda_rational(result.lower_bound()) <= z &&
		    leda_rational(result.upper_bound()) >= z );
    CGAL_FPU_set_rounding_to_infinity();
#endif
    return result;
}

#endif	 // CGAL_IA_LEDA_RATIONAL_H
