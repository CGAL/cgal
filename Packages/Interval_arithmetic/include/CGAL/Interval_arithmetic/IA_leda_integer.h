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
// file          : include/CGAL/Interval_arithmetic/IA_leda_integer.h
// revision      : $Revision$
// revision_date : $Date$
// package       : Interval Arithmetic
// author(s)     : Sylvain Pion
// coordinator   : INRIA Sophia-Antipolis (<Mariette.Yvinec@sophia.inria.fr>)
//
// ============================================================================

#ifndef CGAL_IA_LEDA_INTEGER_H
#define CGAL_IA_LEDA_INTEGER_H

CGAL_BEGIN_NAMESPACE

// We choose the lazy approach, which is good enough: we take the double
// approximation, which is guaranted 1 bit error max, and return an interval
// around this value.  To have something more precise would require access to
// LEDA integer's internal representation, which is not possible without
// modifying LEDA.

inline
Interval_nt_advanced
convert_from_to (const Interval_nt_advanced&, const leda_integer & z)
{
    CGAL_expensive_assertion(FPU_empiric_test() == CGAL_FE_UPWARD);
    FPU_set_cw(CGAL_FE_TONEAREST);
    double approx = CGAL::to_double(z);
    FPU_set_cw(CGAL_FE_UPWARD);
    Interval_nt_advanced result = Interval_nt_advanced(approx) +
	    Interval_nt_advanced::Smallest;
    CGAL_expensive_assertion_code(FPU_set_cw(CGAL_FE_TONEAREST);)
    CGAL_expensive_assertion( leda_integer(result.inf()) <= z &&
		              leda_integer(result.sup()) >= z);
    CGAL_expensive_assertion_code(FPU_set_cw(CGAL_FE_UPWARD);)
    return result;
}

template <>
struct converter<Interval_nt_advanced,leda_integer>
{
    static inline Interval_nt_advanced do_it (const leda_integer & z)
    {
	return convert_from_to(Interval_nt_advanced(), z);
    }
};

CGAL_END_NAMESPACE

#endif // CGAL_IA_LEDA_INTEGER_H
