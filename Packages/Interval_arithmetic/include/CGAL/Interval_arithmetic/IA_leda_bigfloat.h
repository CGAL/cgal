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
// approximation, which is guaranted 1 bit error max, and return an interval
// around this value.  A check for underflow is added.

inline CGAL_Interval_nt_advanced CGAL_to_interval_nt(const leda_bigfloat &z)
{
    const double two_52 = 1.0/(1024.0*1024.0*1024.0*1024.0*1024.0*4.0); //2^-52
    const double approx = CGAL_to_double(z);
    CGAL_Interval_nt_advanced res_ia;

    CGAL_FPU_set_rounding_to_infinity();
    if ((z != 0) && (approx == 0))
    { // We compute the smallest interval, strictly containing zero.
	res_ia = CGAL_Interval_nt_advanced(-two_52,two_52);
	res_ia *= res_ia; // ~ 2^-104
	res_ia *= res_ia; // ~ 2^-208
	res_ia *= res_ia; // ~ 2^-416
	res_ia *= res_ia; // ~ 2^-832
	res_ia *= res_ia; // < 2^-1024
    } else {
	res_ia = CGAL_Interval_nt_advanced (1-two_52,1+two_52)
		 * CGAL_Interval_nt_advanced (approx);
    };
    CGAL_FPU_set_rounding_to_nearest();
    return res_ia;
}

#endif	 // CGAL_IA_LEDA_BIGFLOAT_H
