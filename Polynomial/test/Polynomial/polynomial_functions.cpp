// ============================================================================
//
// Copyright (c) 2001-2008 Max-Planck-Institut Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of EXACUS (http://www.mpi-inf.mpg.de/projects/EXACUS/);
// you may redistribute it under the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with EXACUS.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// ----------------------------------------------------------------------------
//
// Library       : CGAL
// File          : test/Polynomial/polynomial_functions.C
// CGAL_release   : $Name:  $
// Revision      : $Revision$
// Revision_date : $Date$
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

    // coefficient
    assert(CGAL::coefficient(f,0,0) == NT(-1));
    assert(CGAL::coefficient(f,1,0) == NT( 2));
    assert(CGAL::coefficient(f,2,1) == NT( 3));
    assert(CGAL::coefficient(f,1,2) == NT( 1));
    
    // not "existing"
    assert(CGAL::coefficient(f,0,1) == NT( 0));
    assert(CGAL::coefficient(f,2,2) == NT( 0));
    assert(CGAL::coefficient(f,0,3) == NT( 0));
    
    BIPOLY g = BIPOLY(
            UNPOLY(NT(63), NT(50), NT(-37), NT(-85)),
            UNPOLY(NT(49), NT(97), NT(-55)),
            UNPOLY(NT(56), NT(-35)),
            UNPOLY(NT(79))
    );
      
    assert(CGAL::coefficient(g,0,0) == NT( 63));
    assert(CGAL::coefficient(g,1,0) == NT( 50));
    assert(CGAL::coefficient(g,2,0) == NT(-37));
    assert(CGAL::coefficient(g,3,0) == NT(-85));
    assert(CGAL::coefficient(g,0,1) == NT( 49));
    assert(CGAL::coefficient(g,1,1) == NT( 97));
    assert(CGAL::coefficient(g,2,1) == NT(-55));
    assert(CGAL::coefficient(g,0,2) == NT( 56));
    assert(CGAL::coefficient(g,1,2) == NT(-35));
    assert(CGAL::coefficient(g,0,3) == NT( 79));
    // not "existing"
    assert(CGAL::coefficient(f,2,2) == NT( 0));
    assert(CGAL::coefficient(f,0,4) == NT( 0));
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

    // coefficient
    assert(CGAL::coefficient(f1,0,0,0) == NT(-6));
    assert(CGAL::coefficient(f1,0,0,1) == NT( 8));
    assert(CGAL::coefficient(f1,0,1,0) == NT(-3));
    assert(CGAL::coefficient(f1,1,0,0) == NT( 5));
    assert(CGAL::coefficient(f1,0,1,1) == NT( 4));
    assert(CGAL::coefficient(f1,2,0,1) == NT(-2));
    assert(CGAL::coefficient(f1,1,0,2) == NT( 1));
    
    // not "existing"
    assert(CGAL::coefficient(f1,1,0,1) == NT( 0));
    assert(CGAL::coefficient(f1,0,1,2) == NT( 0));
    assert(CGAL::coefficient(f1,1,2,0) == NT( 0));
    assert(CGAL::coefficient(f1,0,0,3) == NT( 0));

    // coefficient
    assert(CGAL::coefficient(f2,0,0,0) == NT(-6));
    assert(CGAL::coefficient(f2,0,0,1) == NT( 8));
    assert(CGAL::coefficient(f2,1,0,0) == NT( 5));
    assert(CGAL::coefficient(f2,2,0,1) == NT(-2));
    assert(CGAL::coefficient(f2,0,1,1) == NT(-5));
    assert(CGAL::coefficient(f2,0,0,3) == NT( 9));
    
    assert(CGAL::coefficient(f2,0,1,0) == NT( 0));
    assert(CGAL::coefficient(f2,2,1,0) == NT( 0));
    assert(CGAL::coefficient(f2,1,1,0) == NT( 0));
    assert(CGAL::coefficient(f2,1,0,1) == NT( 0));
    assert(CGAL::coefficient(f2,1,2,0) == NT( 0));
    assert(CGAL::coefficient(f2,1,1,1) == NT( 0));
    
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
