// ======================================================================
//
// Copyright (c) 1999 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
// 
// release       :
// release_date  :
// 
// file          : src/Interval_arithmetic.C
// revision      : $Revision$
// revision_date : $Date$
// package       : Interval Arithmetic
// author(s)     : Sylvain Pion
//
// coordinator   : INRIA Sophia-Antipolis (<Mariette.Yvinec@sophia.inria.fr>)
//
// ======================================================================
 
#include <CGAL/Interval_arithmetic.h>

CGAL_BEGIN_NAMESPACE

unsigned Interval_nt_advanced::number_of_failures=0;

std::ostream &
operator<< (std::ostream & os, const Interval_nt_advanced & I)
{
    return os << "[" << I.inf() << ";" << I.sup() << "]";
}

std::istream &
operator>> (std::istream & is, Interval_nt_advanced & I)
{
    double d;
    is >> d;
    I = d;
    return is;
}

void force_ieee_double_precision()
{
#ifdef __i386__
    FPU_set_cw(FPU_cw_near);
#endif // __i386__
}

CGAL_END_NAMESPACE
