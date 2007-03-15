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

#include <CGAL/basic.h>

#include <CGAL/Polynomial.h>
#include <CGAL/_test_real_root_isolator.h>

#include <CGAL/Algebraic_kernel_d/Descartes.h>
#include <CGAL/Arithmetic_kernel.h>

template <class AT>
void test_descartes(){
    typedef typename AT::Integer Integer;
    typedef typename AT::Rational Rational;
    {
        typedef ::CGAL::Polynomial<Integer> Polynomial;
        typedef ::CGAL::CGALi::Descartes<Polynomial,Rational> Isolator;
        
        // general test of concept RealRootIsolator
        CGAL::CGALi::test_real_root_isolator<Isolator>();
    }{
        typedef ::CGAL::Polynomial<Rational> Polynomial;
        typedef ::CGAL::CGALi::Descartes<Polynomial,Rational> Isolator;
        // general test of concept RealRootIsolator
        CGAL::CGALi::test_real_root_isolator<Isolator>();
    }    
}
    
int main(){
#ifdef CGAL_USE_LEDA  
    test_descartes<CGAL::LEDA_arithmetic_kernel>();
#endif
#ifdef CGAL_USE_CORE
    test_descartes<CGAL::CORE_arithmetic_kernel>();
#endif
    return EXIT_SUCCESS;
}
// EOF
