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
// file          : include/CGAL/Interval_arithmetic/IA_leda_integer.h
// revision      : $Revision$
// revision_date : $Date$
// package       : Interval Arithmetic
// author(s)     : Sylvain Pion <Sylvain.Pion@sophia.inria.fr>
//
// coordinator   : INRIA Sophia-Antipolis (<Mariette.Yvinec@sophia.inria.fr>)
//
// ============================================================================

#ifndef CGAL_IA_LEDA_INTEGER_H
#define CGAL_IA_LEDA_INTEGER_H

// We choose the lazy approach, which is good enough: we take the double
// approximation, which is guaranted 1 bit error max, and return an interval
// around this value.  To have something more precise would require access to
// LEDA integer's internal representation, which is not possible without
// modifying LEDA.

inline CGAL_Interval_nt_advanced CGAL_convert_to<CGAL_Interval_nt_advanced>
	(const leda_integer &z)
{
    return CGAL_Interval_nt_advanced (CGAL_to_double(z)) +
	   CGAL_Interval_nt_advanced::smallest;
}

#endif	 // CGAL_IA_LEDA_INTEGER_H
