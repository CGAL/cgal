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
// file          : include/CGAL/Interval_arithmetic/IA_Gmpz.h
// revision      : $Revision$
// revision_date : $Date$
// package       : Interval Arithmetic
// author(s)     : Sylvain Pion <Sylvain.Pion@sophia.inria.fr>
//
// coordinator   : INRIA Sophia-Antipolis (<Mariette.Yvinec@sophia.inria.fr>)
//
// ============================================================================

#ifndef CGAL_IA_GMPZ_H
#define CGAL_IA_GMPZ_H

CGAL_BEGIN_NAMESPACE

// We choose the lazy approach, which is good enough: we take the double
// approximation, which is guaranted 1 bit error max (when rounding to
// nearest), and return an interval around this value.
// It should be much faster to have a low level function especially designed
// for that using rounding to infinity.

struct converter<Interval_nt_advanced,Gmpz>
{
    static inline Interval_nt_advanced do_it (const Gmpz & z)
    {
#ifdef CGAL_IA_DEBUG
	CGAL_assertion(FPU_get_cw() == FPU_cw_up);
#endif
	FPU_set_cw(FPU_cw_near);
	double approx = to_double(z);
	FPU_set_cw(FPU_cw_up);
	Interval_nt_advanced result = approx + CGAL_IA_SMALLEST;
#ifdef CGAL_IA_DEBUG
	FPU_set_cw(FPU_cw_near);
	CGAL_assertion(	Gmpz(result.lower_bound()) <= z &&
			Gmpz(result.upper_bound()) >= z);
	FPU_set_cw(FPU_cw_up);
#endif
	return result;
    }
};

CGAL_END_NAMESPACE

#endif	 // CGAL_IA_GMPZ_H
