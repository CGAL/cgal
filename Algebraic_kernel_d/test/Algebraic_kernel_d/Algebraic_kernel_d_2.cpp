// TODO: Add licence
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL:$
// $Id: $
//
//
// Author(s)     : Michael Kerber <mkerber@mpi-inf.mpg.de>
//
// ============================================================================

// code coverage test for Algebraic_curve_kernel_2

// #define CGAL_ACK_DEBUG_FLAG 1

//#define CGAL_AK_ENABLE_DEPRECATED_INTERFACE 1

#include <CGAL/Algebraic_kernel_d/flags.h>

#include <CGAL/Arithmetic_kernel.h>

#include <CGAL/Algebraic_kernel_d_2.h>

#include <CGAL/Sqrt_extension.h>

#include <CGAL/_test_algebraic_kernel_2.h>


template< typename Coefficient >
void test_algebraic_kernel() {
    typedef CGAL::Algebraic_kernel_d_2<Coefficient> Algebraic_kernel_d_2;
    Algebraic_kernel_d_2 ak_2;
    CGAL::test_algebraic_kernel_2<Algebraic_kernel_d_2>(ak_2);
}

int main() {


#ifdef CGAL_HAS_LEDA_ARITHMETIC_KERNEL

#if CGAL_ACK_DEBUG_FLAG
    CGAL_ACK_DEBUG_PRINT << "TESTING LEDA" << std::endl;
#endif
    {

      typedef CGAL::LEDA_arithmetic_kernel AK;
      test_algebraic_kernel<AK::Integer>();
      /*
      test_algebraic_kernel<AK::Rational>();
      test_algebraic_kernel<CGAL::Sqrt_extension<AK::Integer,AK::Integer> >();

      test_algebraic_kernel
        <CGAL::Sqrt_extension<AK::Rational,AK::Rational> >();
      */
    }
#else
    std::cerr << "LEDA tests skipped" << std::endl;
#endif


#ifdef CGAL_HAS_CORE_ARITHMETIC_KERNEL
#if CGAL_ACK_DEBUG_FLAG
    CGAL_ACK_DEBUG_PRINT << "TESTING CORE" << std::endl;
#endif
    {
      typedef CGAL::CORE_arithmetic_kernel AK;
      test_algebraic_kernel<AK::Integer>();
      /*
      test_algebraic_kernel<AK::Rational>();
      test_algebraic_kernel<CGAL::Sqrt_extension<AK::Integer,AK::Integer> >();

      test_algebraic_kernel
        <CGAL::Sqrt_extension<AK::Rational,AK::Rational> >();
      */
    }
#else
    std::cerr << "CORE tests skipped" << std::endl;
#endif


#ifdef CGAL_HAS_GMP_ARITHMETIC_KERNEL
#if CGAL_ACK_DEBUG_FLAG
    CGAL_ACK_DEBUG_PRINT << "TESTING GMP" << std::endl;
#endif
    {
      typedef CGAL::GMP_arithmetic_kernel AK;
      test_algebraic_kernel<AK::Integer>();
      /*
      test_algebraic_kernel<AK::Rational>();
      test_algebraic_kernel<CGAL::Sqrt_extension<AK::Integer,AK::Integer> >();

      test_algebraic_kernel
        <CGAL::Sqrt_extension<AK::Rational,AK::Rational> >();
      */
    }
#else
    std::cerr << "GMP tests skipped" << std::endl;
#endif


    return 0;
}
