// Copyright (c) 1999-2004  Utrecht University (The Netherlands),
// ETH Zurich (Switzerland), Freie Universitaet Berlin (Germany),
// INRIA Sophia-Antipolis (France), Martin-Luther-University Halle-Wittenberg
// (Germany), Max-Planck-Institute Saarbruecken (Germany), RISC Linz (Austria),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Sylvain Pion
 
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

// #include <CGAL/Filtered_exact.h> // only for CGAL_IA_NEW_FILTERS
// But VC++ 7.0 would need some macros defined.
#include <CGAL/FPU.h>

CGAL_BEGIN_NAMESPACE

// Static variables:
#ifdef CGAL_IA_NEW_FILTERS
#include <CGAL/Arithmetic_filter/static_infos/dispatch.h>
#endif


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

#ifdef CGAL_CFG_DENORMALS_COMPILE_BUG
// For compilers which bug on denormalized values at compile time.
// We generate CGAL_IA_MIN_DOUBLE at run time.
namespace {
double init_min_double()
{
    double d = 1;
    double e = 1;
    do {
	d = e;
	e = CGAL_IA_FORCE_TO_DOUBLE(e/2);
    } while (e != 0);
    return d;
}
} // anonymous namespace

namespace CGALi {
double minimin = init_min_double();
}
#endif

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

// needed in order that the test suite passes for Intel7
namespace CGALi {

double zero() { return 0; }

}

CGAL_END_NAMESPACE
