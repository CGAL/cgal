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
// file          : include/CGAL/IA/IA_leda_real.h
// revision      : 1.7
// revision_date : 1 July 1998
// package       : Interval Arithmetic
// author(s)     : Sylvain Pion <Sylvain.Pion@sophia.inria.fr>
//
// coordinator   : INRIA Sophia-Antipolis (<Herve.Bronnimann@sophia.inria.fr>)
//
// ============================================================================

#ifndef CGAL_IA_LEDA_REAL_H
#define CGAL_IA_LEDA_REAL_H

inline CGAL_Interval_nt CGAL_to_interval_nt(const leda_real &z)
{
    const double approx = CGAL_to_double(z);
    const double rel_error = z.get_double_error();
    CGAL_FPU_set_rounding_to_infinity();
    CGAL_Interval_nt_advanced res_ia = ( CGAL_Interval_nt_advanced (1)
	 + CGAL_Interval_nt_advanced(-rel_error,rel_error) )
	* CGAL_Interval_nt_advanced(approx);
    CGAL_FPU_set_rounding_to_nearest();
    return res_ia;
}

#endif	 // CGAL_IA_LEDA_REAL_H
