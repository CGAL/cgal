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

/*! \file NiX/Bitstream_descartes.C
 This is the test file for the class NiX::Bitstream_descartes.

*/

#include <CGAL/basic.h>
#include <cassert>

// include these traits here by 'hand', since not in release 3.3
#include <CGAL/Algebraic_extension_traits.h>
#include <CGAL/Scalar_factor_traits.h>

#include <CGAL/Polynomial.h>

#include <CGAL/_test_real_root_isolator.h>
#include <CGAL/_test_bitstream_descartes.h>

#include <CGAL/Algebraic_kernel_d/Bitstream_descartes.h>
#include <CGAL/Algebraic_kernel_d/Algebraic_real_d_1.h>
#include <CGAL/Arithmetic_kernel.h>

#include <CGAL/Algebraic_kernel_d/Bitstream_descartes_rndl_tree_traits.h>
#include <CGAL/Algebraic_kernel_d/Bitstream_coefficient_kernel.h>

#include <CGAL/Sqrt_extension.h> // used in this file 

template <class AT>
void test_descartes(){
    typedef typename AT::Integer Integer;
    typedef typename AT::Rational Rational;

    {
        typedef CGAL::internal::Bitstream_descartes<
            CGAL::internal::Bitstream_descartes_rndl_tree_traits
            <CGAL::internal::Bitstream_coefficient_kernel<Integer> > > Isolator;
        
        // general test of concept RealRootIsolator
        CGAL::internal::test_real_root_isolator<Isolator>();
    }{
        typedef CGAL::internal::Bitstream_descartes<
            CGAL::internal::Bitstream_descartes_rndl_tree_traits
            <CGAL::internal::Bitstream_coefficient_kernel<Integer> > > Isolator;
        // general test of concept RealRootIsolator
        CGAL::internal::test_real_root_isolator<Isolator>();
    }{
        typedef CGAL::Sqrt_extension<Integer,Integer> EXT;
        typedef typename 
            CGAL::Polynomial_type_generator<EXT,1>::Type Polynomial;
        typedef CGAL::internal::Bitstream_descartes<
            CGAL::internal::Bitstream_descartes_rndl_tree_traits
            <CGAL::internal::Bitstream_coefficient_kernel<EXT> > > Isolator;
        // general test of concept RealRootIsolator
        CGAL::internal::test_real_root_isolator<Isolator>();
    
        std::istringstream is("P[8(0,EXT[1263296571491275162619395552058539312049753537208652637440,42968358109221573436642744060744334362576495937343892480,859])(1,EXT[2207556620237983039471566299950573667219187771717363990528,76852322515373647784745416857429135583058957867416403456,859])(2,EXT[1309275777321138279335848837056750819020098551750287419392,39195296448043974486512553808164864989806318662622537728,859])(3,EXT[86302507833822837152267458208050616275030717971310571520,10003400461933730535898849196215973480541410956350136320,859])(4,EXT[491437197926570070047913809040994179733862944058588160,-823318654055010400576035967724449228967294601081409536,859])(5,EXT[65617171248843379260568930361980285279772972904474353664,-2321236429038088490878641998530459657094667195101177856,859])(6,EXT[-31388640426864731218854617935763549592108582309411053568,1147925677098540153039869220704807888848997763457081344,859])(7,EXT[-9044080753104082029116583549917596926203452476780478464,319259245952387286523925425244746929470795371494383616,859])(8,EXT[7527302869236151900084946597004902733830052515530309632,-256375273905623226297550301204997314964265060295540736,859])]");
    
        Polynomial P;
        is >> P ;
        Isolator isolator(P);
        assert(isolator.number_of_real_roots() == 2 );
        
        typedef CGAL::internal::Algebraic_real_d_1<EXT,Rational> Alg_real;
        Alg_real r0(P,isolator.left_bound(0),isolator.right_bound(0));
        Alg_real r1(P,isolator.left_bound(1),isolator.right_bound(1));
        assert(r0 < r1);
        assert(r0 > isolator.left_bound(0));
        assert(r0 < isolator.right_bound(0));
        assert(r1 > isolator.left_bound(1));
        assert(r1 < isolator.right_bound(1));
    
    }
    CGAL::internal::test_bitstream_descartes<AT>();
    
}
    
int main(){
#ifdef CGAL_HAS_LEDA_ARITHMETIC_KERNEL  
    test_descartes<CGAL::LEDA_arithmetic_kernel>();
#endif

#ifdef CGAL_HAS_CORE_ARITHMETIC_KERNEL
    test_descartes<CGAL::CORE_arithmetic_kernel>();
#endif

#ifdef CGAL_HAS_GMP_ARITHMETIC_KERNEL
    test_descartes<CGAL::GMP_arithmetic_kernel>();
#endif
    return EXIT_SUCCESS;
}
// EOF
