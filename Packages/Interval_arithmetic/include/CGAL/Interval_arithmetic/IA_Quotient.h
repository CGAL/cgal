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

#ifndef IA_QUOTIENT_H
#define IA_QUOTIENT_H

CGAL_BEGIN_NAMESPACE

// We don't know anything about the internal RT type, so there is a risk of
// overflow, but we can't do better than the following trivial conversion.

template <class RT>
struct converter<Interval_nt_advanced,Quotient<RT> >
{
    static inline Interval_nt_advanced do_it (const Quotient<RT> & z)
    {
#ifdef CGAL_IA_DEBUG
	CGAL_assertion(FPU_get_rounding_mode() == FPU_PLUS_INFINITY);
#endif
	return  convert_to<Interval_nt_advanced>(z.numerator()) /
		convert_to<Interval_nt_advanced>(z.denominator());
    }
};

CGAL_END_NAMESPACE

#endif	 // IA_QUOTIENT_H
