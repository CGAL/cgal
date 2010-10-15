// Author(s)     : Michael Hemmer <mhemmer@uni-mainz.de>



/*! \file CGAL/Residue.C
  test for number type modul 
*/

#include <CGAL/basic.h>
#include <cassert>
#include <CGAL/Residue.h>
#include <CGAL/Modular_traits.h>
#include <CGAL/Sqrt_extension.h>
#include <CGAL/Polynomial.h>
#include <CGAL/Lazy_exact_nt.h>
//#include <CGAL/MP_Float.h>
 

#ifdef CGAL_USE_LEDA
#include <CGAL/leda_integer.h>
#include <CGAL/leda_rational.h>
#endif // CGAL_USE_LEDA

#ifdef CGAL_USE_CORE
#include <CGAL/CORE_BigInt.h>
#endif // CGAL_USE_CORE

#ifdef CGAL_USE_GMP
#include <CGAL/Gmpz.h>
#endif

#ifdef CGAL_USE_GMPXX
#include <CGAL/mpz_class.h>
#endif // CGAL_USE_GMP



#include <cstdlib>

#include <boost/type_traits.hpp>

template <class TESTT>
void test_modular_traits(){

        typedef CGAL::Residue Residue;
        typedef CGAL::Modular_traits<TESTT> MT;
        typedef typename MT::Residue_type Residue_type;
        typedef typename MT::Modular_image Modular_image;
        typedef typename MT::Modular_image_representative Modular_image_representative;
        typedef typename MT::Is_modularizable Is_modularizable;
        typedef typename MT::NT NT;
        
        assert(
            !(::boost::is_same<CGAL::Null_functor,Modular_image>::value));
        assert(
            !(::boost::is_same<CGAL::Null_functor,Modular_image_representative>::value));
        assert(
            (::boost::is_same<CGAL::Tag_true,Is_modularizable>::value));
        assert(
            (::boost::is_same<TESTT,NT>::value));
        
        Residue::set_current_prime(7);
        Modular_image modular_image;
        assert(modular_image(TESTT(10)+TESTT(10)) == Residue_type(-1)); 
        assert(modular_image(TESTT(2) *TESTT(10)) == Residue_type(-1)); 
        assert(modular_image(TESTT(20)) == Residue_type(-1)); 
        assert(modular_image(TESTT(20)) == Residue_type(6));   
        assert(modular_image(TESTT(21)) == Residue_type(0));   
        assert(modular_image(TESTT(22)) == Residue_type(1));
        assert(modular_image(TESTT(777777722)) == Residue_type(1));

        Modular_image_representative modular_image_representative;
        assert(modular_image_representative(modular_image(TESTT(20)))
            == TESTT(-1)); 
}

int main()
{ 
  // Enforce IEEE double precision and rounding mode to nearest
  CGAL::Protect_FPU_rounding<true> pfr(CGAL_FE_TONEAREST);
  
    test_modular_traits<int>();
   
#ifdef CGAL_USE_LEDA
    test_modular_traits<leda::integer>();
    test_modular_traits<CGAL::Polynomial< leda::integer > >();
    test_modular_traits<CGAL::Lazy_exact_nt< leda::integer > >();
    test_modular_traits<CGAL::Sqrt_extension< leda::integer , leda::integer > >();
#endif
#ifdef CGAL_USE_CORE
    test_modular_traits<CORE::BigInt>();
    test_modular_traits<CGAL::Polynomial< CORE::BigInt > >();
    test_modular_traits<CGAL::Lazy_exact_nt< CORE::BigInt > >();
    test_modular_traits<CGAL::Sqrt_extension< CORE::BigInt , CORE::BigInt > >();
#endif

#ifdef CGAL_USE_GMP
    test_modular_traits<CGAL::Gmpz>();
#endif 

#ifdef CGAL_USE_GMPXX
    test_modular_traits< mpz_class >();
#endif
    
    // test Sqrt_extension
    test_modular_traits<CGAL::Sqrt_extension< int , int > >();
    assert(
        (!CGAL::Modular_traits<CGAL::Sqrt_extension<double,double> >
            ::Is_modularizable::value));

    // test Polynomial 
    test_modular_traits<CGAL::Polynomial< int > >();
    assert(
        !CGAL::Modular_traits<CGAL::Polynomial<double> >
        ::Is_modularizable::value);

    // test_modular_traits<CGAL::MP_Float >();
    
    test_modular_traits< CGAL::Lazy_exact_nt<int> >();
    assert(
        !CGAL::Modular_traits<CGAL::Lazy_exact_nt< double > >
        ::Is_modularizable::value);
    
    
}
