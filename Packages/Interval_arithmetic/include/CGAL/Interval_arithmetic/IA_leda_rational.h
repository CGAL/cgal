// ============================================================================
//
// Copyright (c) 1998,1999,2000 The CGAL Consortium
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
// ============================================================================

#ifndef CGAL_IA_LEDA_RATIONAL_H
#define CGAL_IA_LEDA_RATIONAL_H

CGAL_BEGIN_NAMESPACE

// For this one, I hope that adding 3 ulps will be enough for a safe
// conversion.  Since LEDA types (except real) don't give information on the
// precision of to_double(), we can't do much...

inline
Interval_nt_advanced
convert_from_to (const Interval_nt_advanced&, const leda_rational & z)
{
    CGAL_expensive_assertion(FPU_empiric_test() == CGAL_FE_UPWARD);
    FPU_set_cw(CGAL_FE_TONEAREST);
    double approx = to_double(z);
    FPU_set_cw(CGAL_FE_UPWARD);

    Interval_nt_advanced result = approx + Interval_nt_advanced::Smallest;
    // We play it safe:
    result += Interval_nt_advanced::Smallest;
    result += Interval_nt_advanced::Smallest;
    CGAL_expensive_assertion_code(FPU_set_cw(CGAL_FE_TONEAREST);)
    CGAL_expensive_assertion( leda_rational(result.inf()) <= z &&
		              leda_rational(result.sup()) >= z );
    CGAL_expensive_assertion_code(FPU_set_cw(CGAL_FE_UPWARD);)
    return result;
}

template <>
struct converter<Interval_nt_advanced,leda_rational>
{
    static inline Interval_nt_advanced do_it (const leda_rational & z)
    {
	return convert_from_to(Interval_nt_advanced(), z);
    }
};

CGAL_END_NAMESPACE

#endif // CGAL_IA_LEDA_RATIONAL_H
