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

template <>
struct converter
{
    static inline Interval_nt_advanced do_it (const leda_rational & z)
    {
#ifdef CGAL_IA_DEBUG
    CGAL_assertion(FPU_get_rounding_mode() == FPU_PLUS_INFINITY);
#endif
    FPU_set_rounding_to_nearest();
    double approx = to_double(z);
    FPU_set_rounding_to_infinity();

    Interval_nt_advanced result = approx + Interval_nt_advanced::smallest;
    // We play it safe:
    result += Interval_nt_advanced::smallest;
    result += Interval_nt_advanced::smallest;
#ifdef CGAL_IA_DEBUG
    FPU_set_rounding_to_nearest();
    CGAL_assertion( leda_rational(result.lower_bound()) <= z &&
		    leda_rational(result.upper_bound()) >= z );
    FPU_set_rounding_to_infinity();
#endif
    return result;
    }
};

CGAL_END_NAMESPACE

#endif	 // CGAL_IA_LEDA_RATIONAL_H
