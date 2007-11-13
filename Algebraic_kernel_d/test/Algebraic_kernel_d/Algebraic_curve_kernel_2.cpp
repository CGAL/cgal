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

#include <CGAL/Benchmark/Benchmark.hpp>
#include <CGAL/Benchmark/Option_parser.hpp>

// reqiured for Kernel_2::decompose tests
#define AcX_CHECK_POLYNOMIALS_FOR_COPRIMABILITY

#include <NiX/Arithmetic_traits.h>
#include <NiX/NT_traits.h>
#include <AcX/Algebraic_curve_2.h>
#include <AcX/Algebraic_curve_pair_2.h>

#include <CGAL/Algebraic_curve_kernel_2.h>
#include <CGAL/_test_algebraic_curve_kernel_2.h>

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
