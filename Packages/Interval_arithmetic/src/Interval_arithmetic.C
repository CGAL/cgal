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

// M$ VC++ doesn't like them yet.
#ifdef CGAL_IA_NEW_FILTERS
#include <CGAL/predicates/kernel_ftC2.h>
#include <CGAL/predicates/kernel_ftC3.h>
#include <CGAL/predicates/sign_of_determinant.h>
#include <CGAL/predicates/Regular_triangulation_ftC2.h>
#include <CGAL/predicates/Regular_triangulation_ftC3.h>
#include <CGAL/predicates/Regular_triangulation_rtH2.h>
#include <CGAL/predicates/Regular_triangulation_rtH3.h>
#endif

#include <CGAL/Arithmetic_filter.h>

CGAL_BEGIN_NAMESPACE

// Static variables:
#ifdef CGAL_IA_NEW_FILTERS
#include <CGAL/Arithmetic_filter/static_infos/dispatch.h>
#endif

unsigned Interval_nt_advanced::number_of_failures = 0;
const Interval_nt_advanced Interval_nt_advanced::Largest (-HUGE_VAL, HUGE_VAL);
const Interval_nt_advanced Interval_nt_advanced::Smallest
             (-CGAL_IA_MIN_DOUBLE, CGAL_IA_MIN_DOUBLE);

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
#if defined __i386__ || defined _MSC_VER || defined __BORLANDC__
    FPU_set_cw(FPU_cw_near);
#endif
}

#ifdef __BORLANDC__
// Borland doesn't initialize the FPU exception mask correctly
// => FP exceptions.
struct Borland_workaround
{
    Borland_workaround() { FPU_set_cw(FPU_cw_near); }
};

static Borland_workaround Borland_workaround_object;
#endif // __BORLANDC__


// The results of 1-epsilon and -1+epsilon are enough
// to detect exactly the current rounding mode.
//                          1-MIN_DOUBLE
//                        +------+-------+
//                        |  1   | 1-ulp |
//               +--------+------+-------+
// -1+MIN_DOUBLE | -1     | near |  -inf |
//               | -1+ulp | +inf |  zero |
//               +--------+------+-------+

FPU_CW_t
FPU_empiric_test()
{
    // If not marked "volatile", the result is false when optimizing
    // because the constants are propagated at compile time !!!
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
