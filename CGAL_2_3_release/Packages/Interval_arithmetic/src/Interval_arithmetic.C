// ======================================================================
//
// Copyright (c) 1999,2000 The CGAL Consortium
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
// coordinator   : INRIA Sophia-Antipolis (<Mariette.Yvinec@sophia.inria.fr>)
//
// ======================================================================
 
#include <CGAL/basic.h>
#include <CGAL/Interval_base.h>

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

#include <CGAL/Filtered_exact.h>

CGAL_BEGIN_NAMESPACE

// Static variables:
#ifdef CGAL_IA_NEW_FILTERS
#include <CGAL/Arithmetic_filter/static_infos/dispatch.h>
#endif

unsigned Interval_base::number_of_failures;

const Interval_base Interval_base::Largest (-HUGE_VAL, HUGE_VAL);
const Interval_base Interval_base::Smallest (-CGAL_IA_MIN_DOUBLE,
                                              CGAL_IA_MIN_DOUBLE);

std::ostream &
operator<< (std::ostream & os, const Interval_base & I)
{
    return os << "[" << I.inf() << ";" << I.sup() << "]";
}

std::istream &
operator>> (std::istream & is, Interval_base & I)
{
    double d;
    is >> d;
    I = d;
    return is;
}

void force_ieee_double_precision()
{
#if defined __i386__ || defined _MSC_VER || defined __BORLANDC__
    FPU_set_cw(CGAL_FE_TONEAREST);
#endif
}

#ifdef __BORLANDC__
// Borland doesn't initialize the FPU exception mask correctly
// => FP exceptions.
struct Borland_workaround
{
    Borland_workaround() { FPU_set_cw(CGAL_FE_TONEAREST); }
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

// I use a global variable here to avoid constant propagation.
double IA_min_double = CGAL_IA_MIN_DOUBLE;

FPU_CW_t
FPU_empiric_test()
{
    double y = 1.0, z = -1.0;
    double ye, ze;
    ye = y - IA_min_double;
    ze = z + IA_min_double;
    if (y == ye && z == ze) return CGAL_FE_TONEAREST;
    if (y == ye) return CGAL_FE_UPWARD;
    if (z == ze) return CGAL_FE_DOWNWARD;
    return CGAL_FE_TOWARDZERO;
}

CGAL_END_NAMESPACE
