// ============================================================================
//
// Copyright (c) 1999 The CGAL Consortium
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
// file          : include/CGAL/Interval_arithmetic/IA_Quotient.h
// revision      : $Revision$
// revision_date : $Date$
// package       : Interval Arithmetic
// author(s)     : Sylvain Pion <Sylvain.Pion@sophia.inria.fr>
//
// coordinator   : INRIA Sophia-Antipolis (<Mariette.Yvinec@sophia.inria.fr>)
//
// ============================================================================

#ifndef CGAL_IA_QUOTIENT_H
#define CGAL_IA_QUOTIENT_H

// We don't know anything about the internal RT type, so there is a risk of
// overflow, but we can't do better than the following trivial conversion.

template <class RT>
inline
CGAL_Interval_nt_advanced
CGAL_convert_to (const CGAL_Quotient<RT> &z, const CGAL_Interval_nt_advanced &)
{
#ifdef CGAL_IA_DEBUG
    CGAL_assertion(CGAL_FPU_get_rounding_mode() == CGAL_FPU_PLUS_INFINITY);
#endif
    return CGAL_convert_to(z.numerator(),   CGAL_Interval_nt_advanced(0)) /
	   CGAL_convert_to(z.denominator(), CGAL_Interval_nt_advanced(0));
}

#endif	 // CGAL_IA_QUOTIENT_H
