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
// revision      : 1.7
// revision_date : 1 July 1998
// package       : Interval Arithmetic
// author(s)     : Sylvain Pion <Sylvain.Pion@sophia.inria.fr>
//
// coordinator   : INRIA Sophia-Antipolis (<Mariette.Yvinec@sophia.inria.fr>)
//
// ============================================================================

#ifndef CGAL_IA_LEDA_RATIONAL_H
#define CGAL_IA_LEDA_RATIONAL_H

// For this one, I prefer not relying on the to_double() member function, as
// it doesn't give any warranty on the precision.

inline CGAL_Interval_nt CGAL_to_interval_nt(const leda_rational &z)
{
    const CGAL_Interval_nt num = CGAL_to_interval_nt(z.numerator());
    const CGAL_Interval_nt den = CGAL_to_interval_nt(z.denominator());
    CGAL_FPU_set_rounding_to_infinity();
    CGAL_Interval_nt_advanced res_ia = num / den;
    CGAL_FPU_set_rounding_to_nearest();
    return res_ia;
}

#endif	 // CGAL_IA_LEDA_RATIONAL_H
