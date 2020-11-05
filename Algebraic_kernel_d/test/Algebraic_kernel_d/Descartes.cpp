// TODO: Add licence
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL:$
// $Id: $
//
//
// Author(s)     :
//
// ============================================================================

// TODO: The comments are all original EXACUS comments and aren't adapted. So
//         they may be wrong now.

/*! \file NiX/Descartes.C
 This is the test file for the class NiX::Descartes.
*/


// include these traits here by 'hand', since not in release 3.3
#include <CGAL/Algebraic_extension_traits.h>
#include <CGAL/Scalar_factor_traits.h>

#include <CGAL/Polynomial_type_generator.h>
#include <CGAL/_test_real_root_isolator.h>

#include <CGAL/Algebraic_kernel_d/Descartes.h>
#include <CGAL/Arithmetic_kernel.h>

template <class AT>
void test_descartes(){
    typedef typename AT::Integer Integer;
    typedef typename AT::Rational Rational;
    {
        typedef typename CGAL::Polynomial_type_generator<Integer,1>::Type
            Polynomial;
        typedef ::CGAL::internal::Descartes<Polynomial,Rational> Isolator;

        // general test of concept RealRootIsolator
        CGAL::internal::test_real_root_isolator<Isolator>();
    }{
        typedef typename CGAL::Polynomial_type_generator<Rational,1>::Type
            Polynomial;
        typedef ::CGAL::internal::Descartes<Polynomial,Rational> Isolator;
        // general test of concept RealRootIsolator
        CGAL::internal::test_real_root_isolator<Isolator>();
    }
}

int main(){
#ifdef CGAL_HAS_LEDA_ARITHMETIC_KERNEL
  std::cout << " TEST AK1 USING LEDA " << std::endl;
  test_descartes< CGAL::LEDA_arithmetic_kernel >();
#endif
#ifdef CGAL_HAS_CORE_ARITHMETIC_KERNEL
  std::cout << " TEST AK1 USING CORE " << std::endl;
  test_descartes< CGAL::CORE_arithmetic_kernel >();
#endif
#ifdef CGAL_HAS_GMP_ARITHMETIC_KERNEL
  std::cout << " TEST AK1 USING GMP " << std::endl;
  test_descartes< CGAL::GMP_arithmetic_kernel >();
#endif

    return EXIT_SUCCESS;
}
// EOF
