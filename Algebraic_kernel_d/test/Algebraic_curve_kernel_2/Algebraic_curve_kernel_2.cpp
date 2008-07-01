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

#ifndef CGAL_ACK_USE_EXACUS
#define CGAL_ACK_USE_EXACUS 0
#endif

#ifndef CGAL_ACK_DEBUG_FLAG
#define CGAL_ACK_DEBUG_FLAG 0
#endif

#include <CGAL/basic.h>

#ifndef CGAL_ACK_2_NO_ALG_REAL_TRAITS_FOR_XY_COORDINATE
#define CGAL_ACK_2_NO_ALG_REAL_TRAITS_FOR_XY_COORDINATE 0
#endif

#include <CGAL/Benchmark/Benchmark.hpp>
#include <CGAL/Benchmark/Option_parser.hpp>

// This is required for the test
#define CGAL_ACK_DONT_CHECK_POLYNOMIALS_FOR_COPRIMALITY 0

#include <CGAL/Arithmetic_kernel.h>
#if CGAL_ACK_USE_EXACUS
#include <AcX/Algebraic_curve_2.h>
#include <AcX/Algebraic_curve_pair_2.h>
#endif

#include <CGAL/Algebraic_kernel_1.h>
#include <CGAL/Algebraic_kernel_d/Algebraic_real_quadratic_refinement_rep_bfi.h>
#include <CGAL/Algebraic_kernel_d/Bitstream_descartes.h>
#include <CGAL/Algebraic_curve_kernel_2.h>

#include <CGAL/Filtered_algebraic_curve_kernel_2.h>

#include <CGAL/_test_algebraic_curve_kernel_2.h>

#include <CGAL/Sqrt_extension.h>

template< class ArithmeticTraits >
void test_algebraic_curve_kernel_2() {

    typedef ArithmeticTraits AT;
    typedef typename AT::Integer Coefficient;
    typedef typename AT::Rational Rational;
      
    typedef CGAL::CGALi::Algebraic_real_quadratic_refinement_rep_bfi
        < Coefficient, Rational > Rep_class;
    typedef CGAL::CGALi::Bitstream_descartes< CGAL::Polynomial< Coefficient >, 
        Rational > Isolator;
    
    typedef CGAL::Algebraic_kernel_1<Coefficient,Rational,Rep_class, Isolator> 
        Algebraic_kernel_1;

#if CGAL_ACK_USE_EXACUS
    typedef AcX::Algebraic_curve_2<Algebraic_kernel_1> Algebraic_curve_2;
    typedef AcX::Algebraic_curve_pair_2<Algebraic_curve_2> 
        Algebraic_curve_pair_2;
    typedef CGAL::Algebraic_curve_kernel_2<Algebraic_curve_pair_2,
        Algebraic_kernel_1> 
        Algebraic_kernel_2;
#else    
    typedef CGAL::Algebraic_curve_kernel_2<Algebraic_kernel_1> 
        Algebraic_kernel_2;
#endif
        
    //std::cout << "Non-filtered kernel..." << std::endl;
    CGAL::CGALi::test_algebraic_curve_kernel_2<Algebraic_kernel_2>();


    //std::cout << "Filtered kernel..." << std::endl;
#if CGAL_ACK_USE_EXACUS
    
    typedef CGAL::Filtered_algebraic_curve_kernel_2<Algebraic_curve_pair_2, 
                                                    Algebraic_kernel_1>
           Filtered_kernel_2;
#else
    typedef CGAL::Filtered_algebraic_curve_kernel_2<Algebraic_kernel_1>
           Filtered_kernel_2;
#endif

    CGAL::CGALi::test_algebraic_curve_kernel_2<Filtered_kernel_2>();
    
}

int main() {

#ifdef CGAL_USE_CORE
    test_algebraic_curve_kernel_2<CGAL::CORE_arithmetic_kernel>();
#else
    std::cerr << "CORE tests skipped" << std::endl;
#endif
#ifdef CGAL_USE_LEDA
    test_algebraic_curve_kernel_2<CGAL::LEDA_arithmetic_kernel>();
#else
    std::cerr << "LEDA tests skipped" << std::endl;
#endif
    
    return 0;
}
