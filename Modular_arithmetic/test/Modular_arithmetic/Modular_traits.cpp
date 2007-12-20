// Author(s)     : Michael Hemmer <mhemmer@uni-mainz.de>

/*! \file CGAL/Modular.C
  test for number type modul 
*/

#include <CGAL/basic.h>
#include <CGAL/Testsuite/assert.h>
#include <CGAL/Modular_traits.h>
#include <CGAL/Sqrt_extension.h>

#ifdef CGAL_USE_LEDA
#include <CGAL/leda_integer.h>
#include <CGAL/leda_rational.h>
#endif // CGAL_USE_LEDA
#ifdef CGAL_USE_CORE
#include <CGAL/CORE_BigInt.h>
#include <CGAL/CORE_BigRat.h>
#endif // CGAL_USE_CORE

#include <cstdlib>

#include <boost/type_traits.hpp>


template <class TESTT>
void test_modular_traits(){

        typedef CGAL::Modular Modular;
        typedef CGAL::Modular_traits<TESTT> MT;
        typedef typename MT::Modular_NT Modular_NT;
        typedef typename MT::Modular_image Modular_image;
        typedef typename MT::Is_modularizable Is_modularizable;
        typedef typename MT::NT NT;
        
        CGAL_test_assert(
                !(::boost::is_same<CGAL::Null_functor,Modular_image>::value));
        CGAL_test_assert(
                (::boost::is_same<CGAL::Tag_true,Is_modularizable>::value));
        CGAL_test_assert(
                (::boost::is_same<TESTT,NT>::value));
        
        Modular::set_current_prime(7);
        Modular_image modular_image;
        CGAL_test_assert(modular_image(TESTT(21)) == Modular_NT(0));   
        CGAL_test_assert(modular_image(TESTT(22)) == Modular_NT(1));
        CGAL_test_assert(modular_image(TESTT(777777722)) == Modular_NT(1));
}

int main()
{   
    test_modular_traits<int>();
#ifdef CGAL_USE_LEDA
    test_modular_traits<leda::integer>();
    test_modular_traits<leda::rational>();
    test_modular_traits<CGAL::Sqrt_extension< leda::integer, leda::integer > >();
#endif
#ifdef CGAL_USE_CORE
    test_modular_traits<CORE::BigInt>();
    test_modular_traits<CORE::BigRat>();
    test_modular_traits<CGAL::Sqrt_extension< CORE::BigInt, CORE::BigInt > >();
#endif

}
