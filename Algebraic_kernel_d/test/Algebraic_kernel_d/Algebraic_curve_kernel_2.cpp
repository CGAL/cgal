// TODO: Add licence
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL:$
// $Id: $
// 
//
// Author(s)     : Pavel Emeliyanenko <asm@mpi-sb.mpg.de>
//
// ============================================================================

// code coverage test for Algebraic_curve_kernel_2

#include <CGAL/basic.h>

#ifndef CGAL_ACK_2_NO_ALG_REAL_TRAITS_FOR_XY_COORDINATE
#define CGAL_ACK_2_NO_ALG_REAL_TRAITS_FOR_XY_COORDINATE 0
#endif

#define NiX_USE_QUADRATIC_REFINEMENT 1
#define NiX_USE_INTERNAL_MODULAR_GCD 1

#include <CGAL/Benchmark/Benchmark.hpp>
#include <CGAL/Benchmark/Option_parser.hpp>

// required for Kernel_2::decompose tests
#define AcX_CHECK_POLYNOMIALS_FOR_COPRIMABILITY 1

#include <NiX/Arithmetic_traits.h>
#include <NiX/NT_traits.h>
#include <AcX/Algebraic_curve_2.h>
#include <AcX/Algebraic_curve_pair_2.h>

#include <CGAL/Algebraic_curve_kernel_2.h>
#include <CGAL/_test_algebraic_curve_kernel_2.h>

#include <CGAL/Sqrt_extension.h>

//#include <CGAL/Arithmetic_kernel.h>

template< class ArithmeticTraits >
void test_algebraic_curve_kernel_2() {

    typedef ArithmeticTraits AT;
    typedef typename AT::Integer Coefficient;
        
    typedef AcX::Algebraic_curve_2<AT> Curve_2;
    typedef AcX::Algebraic_curve_pair_2<Curve_2> Curve_pair_2;
            
    typedef CGAL::Algebraic_kernel_1<Coefficient> Kernel_1;
    typedef CGAL::Algebraic_curve_kernel_2<Curve_pair_2, Kernel_1>
           Kernel_2;

    CGAL::CGALi::test_algebraic_curve_kernel_2<Kernel_2>();
}

int main() {

#ifdef LiS_HAVE_CORE
        typedef NiX::CORE_arithmetic_traits AT;
#else
#ifdef CGAL_USE_LEDA
        typedef NiX::LEDA_arithmetic_traits AT;
#endif
#endif

    test_algebraic_curve_kernel_2<AT>();
    
    return 0;
}
