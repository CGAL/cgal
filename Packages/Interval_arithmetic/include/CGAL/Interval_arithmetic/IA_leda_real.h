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

struct converter<Interval_nt_advanced,leda_real>
{
    static inline Interval_nt_advanced do_it (const leda_real & z)
    {
#ifdef CGAL_IA_DEBUG
    CGAL_assertion(FPU_get_cw() == FPU_cw_up);
#endif
    FPU_set_cw(FPU_cw_near);
    double approx = to_double(z);
    double rel_error = z.get_double_error();
    FPU_set_cw(FPU_cw_up);
    Interval_nt_advanced result = approx
	* ( Interval_nt_advanced(-rel_error,rel_error) + 1 );
#ifdef CGAL_IA_DEBUG
    FPU_set_cw(FPU_cw_near);
    CGAL_assertion( leda_real(result.lower_bound()) <= z &&
		    leda_real(result.upper_bound()) >= z );
    FPU_set_cw(FPU_cw_up);
#endif
    return result;
    }
};

CGAL_END_NAMESPACE

#endif	 // CGAL_IA_LEDA_REAL_H
