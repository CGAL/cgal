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
 
#include <CGAL/basic.h>

#include <CGAL/predicates_on_ftC2.h>
#include <CGAL/predicates_on_ftC3.h>
#include <CGAL/predicates/sign_of_determinant.h>
#include <CGAL/predicates/Regular_triangulation_ftC2.h>
#include <CGAL/predicates/Regular_triangulation_ftC3.h>
#include <CGAL/predicates/Regular_triangulation_rtH2.h>
#include <CGAL/predicates/Regular_triangulation_rtH3.h>

#include <CGAL/Arithmetic_filter.h>

CGAL_BEGIN_NAMESPACE

// Static variables:
#include <CGAL/Arithmetic_filter/static_infos/dispatch.h>

unsigned Interval_nt_advanced::number_of_failures = 0;
bool     Interval_nt_advanced::want_exceptions    = true;

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

// The results of 1-epsilon and -1+epsilon are enough
// to detect exactly the rounding mode.
// ----------------------------------------------------
// rounding mode:        +inf    -inf    0       nearest
// ----------------------------------------------------
//  1-MIN_DOUBLE         1       1-ulp   1-ulp   1
// -1+MIN_DOUBLE        -1+ulp  -1      -1+ulp  -1
// ----------------------------------------------------

FPU_CW_t FPU_empiric_test()
{
    // If not marked "volatile", the result is false when optimizing
    // because the constants are pre-computed at compile time !!!
    volatile const double m = CGAL_IA_MIN_DOUBLE;
    const double y = 1.0, z = -1.0;
    double ye, ze;
    ye = y - m;
    ze = z + m;
    if (y == ye && z == ze) return FPU_cw_near;
    if (y == ye) return FPU_cw_up;
    if (z == ze) return FPU_cw_down;
    return FPU_cw_zero;
}

CGAL_END_NAMESPACE
