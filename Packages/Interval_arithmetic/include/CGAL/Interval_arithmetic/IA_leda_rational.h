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

CGAL_BEGIN_NAMESPACE

// For this one, I hope that adding 3 ulps will be enough for a safe
// conversion.  Since LEDA types (except real) don't give information on the
// precision of to_double(), we can't do much...

struct converter<Interval_nt_advanced,leda_rational>
{
    static inline Interval_nt_advanced do_it (const leda_rational & z)
    {
#ifdef CGAL_IA_DEBUG
    CGAL_assertion(FPU_get_cw() == FPU_cw_up);
#endif
    FPU_set_cw(FPU_cw_near);
    double approx = to_double(z);
    FPU_set_cw(FPU_cw_up);

    Interval_nt_advanced result = approx + CGAL_IA_SMALLEST;
    // We play it safe:
    result += CGAL_IA_SMALLEST;
    result += CGAL_IA_SMALLEST;
#ifdef CGAL_IA_DEBUG
    FPU_set_cw(FPU_cw_near);
    CGAL_assertion( leda_rational(result.lower_bound()) <= z &&
		    leda_rational(result.upper_bound()) >= z );
    FPU_set_cw(FPU_cw_up);
#endif
    return result;
    }
};

CGAL_END_NAMESPACE

#endif	 // CGAL_IA_LEDA_RATIONAL_H
