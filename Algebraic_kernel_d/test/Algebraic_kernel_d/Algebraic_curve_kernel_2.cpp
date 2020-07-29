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

#define CGAL_ACK_DEBUG_FLAG 0

#define CGAL_AK_ENABLE_DEPRECATED_INTERFACE 1

#include <CGAL/Algebraic_kernel_d/flags.h>

#include <CGAL/Arithmetic_kernel.h>

#if CGAL_ACK_USE_EXACUS
#include <AcX/Algebraic_curve_2.h>
#include <AcX/Algebraic_curve_pair_2.h>
#endif

#include <CGAL/Algebraic_kernel_d_1.h>
#include <CGAL/Algebraic_kernel_d/Algebraic_real_quadratic_refinement_rep_bfi.h>
#include <CGAL/Algebraic_kernel_d/Bitstream_descartes.h>
#include <CGAL/Algebraic_kernel_d/Algebraic_curve_kernel_2.h>

#if CGAL_ACK_WITH_FILTERED_KERNEL
#include <CGAL/Filtered_algebraic_curve_kernel_2.h>
#endif

#include <CGAL/_test_algebraic_curve_kernel_2.h>

#include <CGAL/Sqrt_extension.h>

template< class ArithmeticTraits >
void test_algebraic_curve_kernel_2() {

    typedef ArithmeticTraits AT;
    typedef typename AT::Integer Coefficient;
    typedef typename AT::Rational Rational;

    typedef CGAL::internal::Algebraic_real_quadratic_refinement_rep_bfi
        < Coefficient, Rational > Rep_class;
    typedef CGAL::internal::Bitstream_descartes<
        CGAL::internal::Bitstream_descartes_rndl_tree_traits<
        CGAL::internal::Bitstream_coefficient_kernel<Coefficient > > >
        Isolator;

    typedef CGAL::Algebraic_kernel_d_1<Coefficient,Rational,Rep_class, Isolator>
        Algebraic_kernel_d_1;

#if CGAL_ACK_USE_EXACUS
    typedef AcX::Algebraic_curve_2<Algebraic_kernel_d_1> Algebraic_curve_2;
    typedef AcX::Algebraic_curve_pair_2<Algebraic_curve_2>
        Algebraic_curve_pair_2;
    typedef CGAL::Algebraic_curve_kernel_2<Algebraic_curve_pair_2,
        Algebraic_kernel_d_1>
        Algebraic_kernel_d_2;
#else
    typedef CGAL::Algebraic_curve_kernel_2<Algebraic_kernel_d_1>
        Algebraic_kernel_d_2;
#endif

    //std::cout << "Non-filtered kernel..." << std::endl;
    CGAL::internal::test_algebraic_curve_kernel_2<Algebraic_kernel_d_2>();

#if CGAL_ACK_WITH_FILTERED_KERNEL
    //std::cout << "Filtered kernel..." << std::endl;
#if CGAL_ACK_USE_EXACUS

    typedef CGAL::Filtered_algebraic_curve_kernel_2<Algebraic_curve_pair_2,
                                                    Algebraic_kernel_d_1>
           Filtered_kernel_2;
#else
    typedef CGAL::Filtered_algebraic_curve_kernel_2<Algebraic_kernel_d_1>
           Filtered_kernel_2;
#endif

    CGAL::internal::test_algebraic_curve_kernel_2<Filtered_kernel_2>();

#endif
}

int main() {


#ifdef CGAL_HAS_LEDA_ARITHMETIC_KERNEL

#if CGAL_ACK_DEBUG_FLAG
    CGAL_ACK_DEBUG_PRINT << "TESTING LEDA" << std::endl;
#endif
    test_algebraic_curve_kernel_2<CGAL::LEDA_arithmetic_kernel>();
#else
    std::cerr << "LEDA tests skipped" << std::endl;
#endif


#ifdef CGAL_HAS_CORE_ARITHMETIC_KERNEL
#if CGAL_ACK_DEBUG_FLAG
    CGAL_ACK_DEBUG_PRINT << "TESTING CORE" << std::endl;
#endif

    test_algebraic_curve_kernel_2<CGAL::CORE_arithmetic_kernel>();
#else
    std::cerr << "CORE tests skipped" << std::endl;
#endif


#ifdef CGAL_HAS_GMP_ARITHMETIC_KERNEL
#if CGAL_ACK_DEBUG_FLAG
    CGAL_ACK_DEBUG_PRINT << "TESTING GMP" << std::endl;
#endif
    test_algebraic_curve_kernel_2<CGAL::GMP_arithmetic_kernel>();
#else
    std::cerr << "GMP tests skipped" << std::endl;
#endif


    return 0;
}
