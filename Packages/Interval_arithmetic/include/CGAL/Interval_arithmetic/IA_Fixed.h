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

CGAL_BEGIN_NAMESPACE

// The Fixed is in fact a float => trivial conversion.

inline
Interval_nt_advanced
convert_from_to (const Interval_nt_advanced&, const Fixed_precision_nt & z)
{
    return to_double(z);
}

#ifndef CGAL_CFG_NO_EXPLICIT_TEMPLATE_FUNCTION_ARGUMENT_SPECIFICATION
template <>
struct converter<Interval_nt_advanced,Fixed_precision_nt>
{
    static inline Interval_nt_advanced do_it (const Fixed_precision_nt & z)
    { return to_double(z); }
};
#endif // CGAL_CFG_NO_EXPLICIT_TEMPLATE_FUNCTION_ARGUMENT_SPECIFICATION


CGAL_END_NAMESPACE

#endif	 // CGAL_IA_FIXED_H
