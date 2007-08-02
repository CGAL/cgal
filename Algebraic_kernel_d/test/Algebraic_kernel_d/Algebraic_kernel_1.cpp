// TODO: Add licence
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL:$
// $Id: $
// 
//
// Author(s)     : Sebastian Limbach <slimbach@mpi-inf.mpg.de>
//                 Michael Hemmer    <hemmer@mpi-inf.mpg.de>
//
// ============================================================================

// Test of Algebraic_kernel

#include <CGAL/basic.h>
#include <CGAL/Algebraic_kernel_1.h>
#include <CGAL/Algebraic_kernel_d/Algebraic_real_rep_bfi.h>
#include <CGAL/_test_algebraic_kernel_1.h>

#include <CGAL/Arithmetic_kernel.h>


template< class Coefficient_, class Boundary_, class RepClass >
void test_algebraic_kernel_intern() {
    typedef Coefficient_ Coefficient;
    typedef Boundary_    Boundary;
    typedef RepClass     Rep_class;

    typedef CGAL::Algebraic_kernel_1< Coefficient, Boundary, Rep_class > 
                                                                Kernel;
    typedef CGAL::CGALi::Algebraic_real_pure< 
             Coefficient, 
             Boundary,
            ::CGAL::Handle_policy_no_union,     
            Rep_class >                                         Algebraic_real_1;
    typedef CGAL::Polynomial< Coefficient >                     Polynomial_1;
    typedef CGAL::CGALi::Descartes< Polynomial_1, Boundary >    Isolator;

    CGAL::CGALi::test_algebraic_kernel_1< 
        Kernel, Algebraic_real_1, Isolator, Coefficient, Polynomial_1, Boundary >();    
}

template< class ArithmeticKernel >
void test_algebraic_kernel() {
    typedef ArithmeticKernel AK;
    
    // Test with integer
    test_algebraic_kernel_intern< typename AK::Integer,
                                  typename AK::Rational,
                                  CGAL::CGALi::Algebraic_real_rep< typename AK::Integer,
                                                                   typename AK::Rational > >();
   // Test with integer and bfi
    test_algebraic_kernel_intern< typename AK::Integer,
                                  typename AK::Rational,
                                  CGAL::CGALi::Algebraic_real_rep_bfi< typename AK::Integer,
                                                                       typename AK::Rational > >();
                                                                       
    // Test with rational
    test_algebraic_kernel_intern< typename AK::Rational,
                                  typename AK::Rational,
                                  CGAL::CGALi::Algebraic_real_rep< typename AK::Rational,
                                                                   typename AK::Rational > >();

    // Test with rational and bfi
    test_algebraic_kernel_intern< typename AK::Rational,
                                  typename AK::Rational,
                                  CGAL::CGALi::Algebraic_real_rep_bfi< typename AK::Rational,
                                                                       typename AK::Rational > >();
                                                                       
    // Test with Sqrt_extension
    test_algebraic_kernel_intern< CGAL::Sqrt_extension< typename AK::Integer, typename AK::Integer>,
                                  typename AK::Rational,
                                  CGAL::CGALi::Algebraic_real_rep< CGAL::Sqrt_extension< typename AK::Integer, typename AK::Integer >,
                                                                   typename AK::Rational > >();

}

int main() {
#ifdef CGAL_USE_LEDA
    test_algebraic_kernel< CGAL::LEDA_arithmetic_kernel >();
#endif
    
    return 0;
}
