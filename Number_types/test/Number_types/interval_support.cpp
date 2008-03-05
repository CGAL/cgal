


#include <CGAL/basic.h>
#include <CGAL/Testsuite/assert.h>
#include <CGAL/interval_support.h>

#ifdef CGAL_USE_CORE
#include <CGAL/core_interval_support.h>
#endif

#ifdef CGAL_USE_LEDA
#include <CGAL/leda_interval_support.h>
#endif


// TODO: separat BFI functors from interval functors. 

template<class BFI>
void generic_test_bigfloat_interval(){
    
    typedef CGAL::Bigfloat_interval_traits<BFI> BFIT;
    typedef typename BFIT::BF BF;
    
    //TODO: move this into an new Interval_traits
    // get_significant_bits
    CGAL_test_assert(CGAL::get_significant_bits(BFI(3)) == 2);
    // upper
    CGAL_test_assert(CGAL::upper(BFI(3)) == BF(3));
    // lower
    CGAL_test_assert(CGAL::lower(BFI(3)) == BF(3));
    
    // TODO: rm BFI() within call 
    // get/set_precsion
    long precision = CGAL::get_precision(BFI());
    CGAL_test_assert(CGAL::set_precision(BFI(),15) == precision);
    CGAL_test_assert(CGAL::set_precision(BFI(),precision) == 15);
    
    // hull
    CGAL_test_assert(CGAL::lower(CGAL::hull(BFI(2),BFI(5))) >= BF(1));
    CGAL_test_assert(CGAL::lower(CGAL::hull(BFI(2),BFI(5))) <= BF(2));
    CGAL_test_assert(CGAL::upper(CGAL::hull(BFI(2),BFI(5))) >= BF(5));
    CGAL_test_assert(CGAL::upper(CGAL::hull(BFI(2),BFI(5))) <= BF(6));

    // in_zero
    CGAL_test_assert(CGAL::in_zero(CGAL::hull(BFI( 2),BFI(3))) == false);
    CGAL_test_assert(CGAL::in_zero(CGAL::hull(BFI(-2),BFI(3))) == true);
    
    // overlap
    BFI hull_2_5 = CGAL::hull(BFI(2),BFI(5));
    CGAL_test_assert(CGAL::overlap(hull_2_5, CGAL::hull(BFI(6),BFI(7))) == false);
    CGAL_test_assert(CGAL::overlap(hull_2_5, CGAL::hull(BFI(5),BFI(6))) == true);
    CGAL_test_assert(CGAL::overlap(hull_2_5, CGAL::hull(BFI(4),BFI(5))) == true);
    CGAL_test_assert(CGAL::overlap(hull_2_5, CGAL::hull(BFI(3),BFI(4))) == true);
    CGAL_test_assert(CGAL::overlap(hull_2_5, CGAL::hull(BFI(2),BFI(3))) == true);
    CGAL_test_assert(CGAL::overlap(hull_2_5, CGAL::hull(BFI(1),BFI(2))) == true);
    CGAL_test_assert(CGAL::overlap(hull_2_5, CGAL::hull(BFI(0),BFI(1))) == false);

    // singelton
    CGAL_test_assert(CGAL::singleton(CGAL::hull(BFI(2),BFI(2))) == true);
    CGAL_test_assert(CGAL::singleton(CGAL::hull(BFI(2),BFI(3))) == false);

    // width
    CGAL_test_assert(CGAL::width(CGAL::hull(BFI(2),BFI(2))) == BF(0));
    CGAL_test_assert(CGAL::width(CGAL::hull(BFI(2),BFI(3))) == BF(1));   

    typedef typename CGAL::Get_arithmetic_kernel<BFI>::Arithmetic_kernel AK;
    typedef typename AK::Integer Integer; 
    typedef typename AK::Rational Rational; 
    typedef typename AK::Field_with_sqrt FWS;
    
    CGAL_test_assert(CGAL::convert_to_bfi(Integer(1)) == BFI(1));
    CGAL_test_assert(CGAL::convert_to_bfi(Rational(1)) == BFI(1));
    CGAL_test_assert(CGAL::convert_to_bfi(FWS(1)) == BFI(1));   
}


int main(){

#ifdef CGAL_USE_CORE
    // generic_test_bigfloat_interval<CORE::BigFloat>();
#endif 

#ifdef CGAL_USE_LEDA
    generic_test_bigfloat_interval<CGAL::leda_bigfloat_interval>();
#endif 

}
