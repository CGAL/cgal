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
// file          : include/CGAL/Interval_arithmetic/IA_leda_integer.h
// revision      : 1.7
// revision_date : 1 July 1998
// package       : Interval Arithmetic
// author(s)     : Sylvain Pion <Sylvain.Pion@sophia.inria.fr>
//
// coordinator   : INRIA Sophia-Antipolis (<Mariette.Yvinec@sophia.inria.fr>)
//
// ============================================================================

#ifndef CGAL_IA_LEDA_INTEGER_H
#define CGAL_IA_LEDA_INTEGER_H

// We choose the lazy approach, which is good enough: we take the double
// approximation, which is guaranted 1 bit error max, and return an interval
// around this value.  To have something more precise would require access to
// LEDA integer's internal representation, which is not possible without
// modifying LEDA.

inline CGAL_Interval_nt CGAL_to_interval_nt(const leda_integer &z)
{
    const double two_52 = 1.0/(1024.0*1024.0*1024.0*1024.0*1024.0*4.0); // 2^-52
    const double approx = CGAL_to_double(z);
    CGAL_FPU_set_rounding_to_infinity();
    CGAL_Interval_nt_advanced res_ia =
	CGAL_Interval_nt_advanced (1-two_52,1+two_52) *
	CGAL_Interval_nt_advanced(approx);
    CGAL_FPU_set_rounding_to_nearest();
    return res_ia;
}

#endif	 // CGAL_IA_LEDA_INTEGER_H
