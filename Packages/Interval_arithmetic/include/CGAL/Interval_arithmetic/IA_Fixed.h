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
// file          : include/CGAL/Interval_arithmetic/IA_Fixed.h
// revision      : $Revision$
// revision_date : $Date$
// package       : Interval Arithmetic
// author(s)     : Sylvain Pion <Sylvain.Pion@sophia.inria.fr>
//
// coordinator   : INRIA Sophia-Antipolis (<Mariette.Yvinec@sophia.inria.fr>)
// ============================================================================

#ifndef CGAL_IA_FIXED_H
#define CGAL_IA_FIXED_H

CGAL_BEGIN_NAMESPACE

// Fixed is in fact a float => trivial conversion.

inline
Interval_nt_advanced
convert_from_to (const Interval_nt_advanced&, const Fixed_precision_nt & z)
{
    return to_double(z);
}

template <>
struct converter<Interval_nt_advanced,Fixed_precision_nt>
{
    static inline Interval_nt_advanced do_it (const Fixed_precision_nt & z)
    { return to_double(z); }
};

CGAL_END_NAMESPACE

#endif // CGAL_IA_FIXED_H
