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
convert_to (const leda_real &z, const Interval_nt_advanced &)
{
#ifdef CGAL_IA_DEBUG
    CGAL_assertion(FPU_get_rounding_mode() == FPU_PLUS_INFINITY);
#endif
    FPU_set_rounding_to_nearest();
    const double approx = to_double(z);
    const double rel_error = z.get_double_error();
    FPU_set_rounding_to_infinity();
    const Interval_nt_advanced result =
	( Interval_nt_advanced(-rel_error,rel_error) + 1 )
	* Interval_nt_advanced(approx);
#ifdef CGAL_IA_DEBUG
    FPU_set_rounding_to_nearest();
    CGAL_assertion( leda_real(result.lower_bound()) <= z &&
		    leda_real(result.upper_bound()) >= z );
    FPU_set_rounding_to_infinity();
#endif
    return result;
}

CGAL_END_NAMESPACE

#endif	 // CGAL_IA_LEDA_REAL_H
