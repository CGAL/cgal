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
// file          : include/CGAL/Interval_arithmetic/IA_Fixed.h
// revision      : $Revision$
// revision_date : $Date$
// package       : Interval Arithmetic
// author(s)     : Sylvain Pion <Sylvain.Pion@sophia.inria.fr>
//
// coordinator   : INRIA Sophia-Antipolis (<Mariette.Yvinec@sophia.inria.fr>)
//
// ============================================================================

#ifndef CGAL_IA_FIXED_H
#define CGAL_IA_FIXED_H

// The Fixed is in fact a float => exact conversion.

inline
CGAL_Interval_nt_advanced
CGAL_convert_to (const CGAL_Fixed_precision_nt &z,
		 const CGAL_Interval_nt_advanced &)
{
    return CGAL_Interval_nt_advanced (z.to_double());
}

#endif	 // CGAL_IA_FIXED_H
