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

// For this one, I prefer not relying on the to_double() member function, as
// it doesn't give any warranty on the precision.

inline CGAL_Interval_nt_advanced CGAL_convert_to<CGAL_Interval_nt_advanced>
	(const leda_rational &z)
{
    return CGAL_convert_to<CGAL_Interval_nt_advanced>(z.numerator())
	/  CGAL_convert_to<CGAL_Interval_nt_advanced>(z.denominator());
}

#endif	 // CGAL_IA_LEDA_RATIONAL_H
