// ============================================================================
//
// Copyright (c) 2001-2008 Max-Planck-Institut Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of EXACUS (http://www.mpi-inf.mpg.de/projects/EXACUS/).
//
// ----------------------------------------------------------------------------
//
// Library       : CGAL
// File          : test/Polynomial/polynomial_functions.C
// CGAL_release   : $Name:  $
// Revision      : $Revision$
// Revision_date : $Date$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Eric Berberich <eric@mpi-inf.mpg.de>
//                 Arno Eigenwillig <arno@mpi-inf.mpg.de>
//
// ============================================================================

/*! \file polynomial_functions.C
  test for functions related to polynomials
*/

#include <CGAL/config.h>
#include <CGAL/Arithmetic_kernel.h>
#include <CGAL/CORE_arithmetic_kernel.h>
#include <CGAL/LEDA_arithmetic_kernel.h>
#include <CGAL/Polynomial.h>
#include <CGAL/Polynomial/polynomial_functions.h>
#include <cassert>
#include <cstdlib>

typedef int NT;
typedef CGAL::Polynomial<NT> UNPOLY;
typedef CGAL::Polynomial<UNPOLY> BIPOLY;
typedef CGAL::Polynomial<BIPOLY> TRIPOLY;

template <class AK>
void bivariate_polynomial_test() {

    BIPOLY f0;

    // f = xy^2 + 3x^2y + 2x - 1
    BIPOLY f = BIPOLY(UNPOLY(NT(-1), NT(2)), UNPOLY(NT(0), NT(0), NT(3)),
                      UNPOLY(NT(0), NT(1))
    );

    // differentiation
    // fx = y^2 + 6xy + 2
    BIPOLY fx = BIPOLY(UNPOLY(NT(2)), UNPOLY(NT(0), NT(6)),
            UNPOLY(NT(1))
    );

    // LCC
    assert(CGAL::check_leadcoeff(f0) == true);
    assert(CGAL::check_leadcoeff(f) == false);
    assert(CGAL::check_leadcoeff(fx) == true);
}


template <class AK>
void trivariate_polynomial_test() {


    TRIPOLY f0;
    assert(CGAL::check_leadcoeff(f0) == true);

    // f1 = xz^2 + 4yz - 2x^2z + 5x - 3y + 8z -6
    TRIPOLY f1 = TRIPOLY(
            BIPOLY(UNPOLY(NT(-6), NT(5)), UNPOLY(NT(-3))),
            BIPOLY(UNPOLY(NT(8), NT(0), NT(-2)), UNPOLY(NT(4))),
            BIPOLY(UNPOLY(NT(0),NT(1)))
    );

    // f2 = 9z^3 + -5yz - 2x^2z + 5x + 8z -6
    TRIPOLY f2 = TRIPOLY(
            BIPOLY(UNPOLY(NT(-6), NT(5))),
            BIPOLY(UNPOLY(NT(8), NT(0), NT(-2)), UNPOLY(NT(-5))),
            BIPOLY(UNPOLY(NT(0))),
            BIPOLY(UNPOLY(NT(9)))
    );

    // LCC
    assert(CGAL::check_leadcoeff(f1) == false);
    assert(CGAL::check_leadcoeff(f2) == true);
}


template <class AK>
void polynomial_functions_test() {
    ::CGAL::set_pretty_mode(std::cout);

    //std::cout<<" univariate "<<std::endl;
    //univariate_polynomial_test<AK>();
    std::cout<<" bivariate "<<std::endl;
    bivariate_polynomial_test<AK>();
    std::cout<<" trivariate "<<std::endl;
    trivariate_polynomial_test<AK>();
}


int main(){
#ifdef CGAL_USE_LEDA
    polynomial_functions_test<CGAL::LEDA_arithmetic_kernel>();
#endif // CGAL_USE_LEDA
#ifdef CGAL_USE_CORE
    polynomial_functions_test<CGAL::CORE_arithmetic_kernel>();
#endif // CGAL_USE_CORE
}
// EOF
