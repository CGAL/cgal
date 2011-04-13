// ============================================================================
//
// Copyright (c) 2001 The CGAL Consortium
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
// file          : include/CGAL/NT_converter.h
// revision      : $Revision$
// revision_date : $Date$
// package       : Number_types
// author(s)     : Sylvain Pion
// coordinator   : INRIA Sophia-Antipolis (<Mariette.Yvinec@sophia.inria.fr>)
//
// ============================================================================

#ifndef CGAL_NT_CONVERTER_H
#define CGAL_NT_CONVERTER_H

CGAL_BEGIN_NAMESPACE

// A number type converter usable as default, using the conversion operator.

template < class NT1, class NT2 >
struct NT_converter
{
    NT2
    operator()(const NT1 &a) const
    {
        return NT2(a);
    }
};


// A number type converter to Interval_base, using to_interval().

template < class NT >
struct Interval_converter
{
    Interval_base
    operator()(const NT &a) const
    {
        return to_interval(a);
    }
};

CGAL_END_NAMESPACE

#endif // CGAL_NT_CONVERTER_H
