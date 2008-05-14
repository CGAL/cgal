#include <CGAL/basic.h>
#include <cassert>
#include <CGAL/Arithmetic_kernel.h>
#include <CGAL/convert_to_bfi.h>

#include <CGAL/Test/_test_interval.h>

template<class Interval_>
void generic_test_interval(){
   
    typedef typename CGAL::Interval_traits<Interval_>::Self IT;
    typedef typename IT::Interval Interval;
    typedef typename IT::Boundary Boundary; 
    
    const typename IT::Construct construct = typename IT::Construct();
    const typename IT::Lower lower = typename IT::Lower();
    const typename IT::Upper upper = typename IT::Upper();
    const typename IT::Width width = typename IT::Width();
    const typename IT::Median median = typename IT::Median();
    const typename IT::Norm norm = typename IT::Norm();

    // const typename IT::Empty empty(); !for CORE
    const typename IT::Singleton singleton = typename IT::Singleton();
    const typename IT::Zero_in zero_in = typename IT::Zero_in();
    const typename IT::In in = typename IT::In();
    const typename IT::Equal equal = typename IT::Equal();
    const typename IT::Overlap overlap = typename IT::Overlap();
    const typename IT::Subset subset = typename IT::Subset();
    const typename IT::Proper_subset proper_subset = typename IT::Proper_subset();
    
    const typename IT::Intersection intersection = typename IT::Intersection();
    const typename IT::Hull hull = typename IT::Hull();
    
    Interval a(construct(Boundary(-7),Boundary(-5)));
    Interval b(construct(Boundary(0),Boundary(4)));
    Interval c(construct(Boundary(2),Boundary(6)));

    assert(lower(a)  == Boundary(-7));
    assert(upper(a)  == Boundary(-5));
    assert(lower(b)  == Boundary( 0));
    assert(upper(b)  == Boundary( 4));
    assert(lower(c)  == Boundary( 2));
    assert(upper(c)  == Boundary( 6));

    assert(width(a)  == Boundary( 2));
    assert(median(a) == Boundary(-6));
    assert(norm(a)   == Boundary( 7));
    
    // assert(!empty(a));
    assert( singleton(Interval(1)));
    assert(!singleton(a));
    assert(!singleton(b));
    assert(!singleton(c));
    
    assert(!zero_in(Interval(1)));
    assert( zero_in(Interval(0)));
    assert(!zero_in(a));
    assert( zero_in(b));
    assert(!zero_in(c));

//########
    // to be remove again
    assert(!CGAL::in_zero(Interval(1)));
    assert( CGAL::in_zero(Interval(0)));
    assert(!CGAL::in_zero(a));
    assert( CGAL::in_zero(b));
    assert(!CGAL::in_zero(c));
//#########
    
    assert(!in(Boundary( 3),a));
    assert( in(Boundary(-7),a));
    
    
    assert( equal(a,a));
    assert( equal(b,b));
    assert( equal(c,c));
    assert(!equal(a,b));
    assert(!equal(a,c));

    assert(!overlap(a,b));
    assert( overlap(b,c));
    Interval I25 = construct(Boundary(2),Boundary(5));
    assert(overlap(I25, construct(Boundary(6),Boundary(7))) == false);
    assert(overlap(I25, construct(Boundary(5),Boundary(6))) == true);
    assert(overlap(I25, construct(Boundary(4),Boundary(5))) == true);
    assert(overlap(I25, construct(Boundary(3),Boundary(4))) == true);
    assert(overlap(I25, construct(Boundary(2),Boundary(3))) == true);
    assert(overlap(I25, construct(Boundary(1),Boundary(2))) == true);
    assert(overlap(I25, construct(Boundary(0),Boundary(1))) == false);
    
    assert(!subset(a,b));
    assert( subset(a,a));
    assert( subset(Interval(-6),a));
    
    assert(!proper_subset(a,b));
    assert(!proper_subset(a,a));
    assert( proper_subset(Interval(-6),a));
    
    // assert( empty(intersection(a,b)));
    assert( lower(intersection(b,c)) == Boundary(2));
    assert( upper(intersection(b,c)) == Boundary(4));
    // this part chages in case we allow empty intersection 
    // which seems to be not possible for CORE::BigFloat as Interval 
    try{
        try{
            intersection(a,b);
            assert(false); // it should not reach this 
        }
        catch(CGAL::Exception_intersection_is_empty){} // it throws the right exception 
    }catch(...){
        assert(false); // seems to be the wrong exception
    }
    
    // hull
    assert(lower(hull(b,c)) == Boundary(0));
    assert(upper(hull(b,c)) == Boundary(6));  
    assert(lower(hull(Interval(2),Interval(5))) >= Boundary(1));
    assert(lower(hull(Interval(2),Interval(5))) <= Boundary(2));
    assert(upper(hull(Interval(2),Interval(5))) >= Boundary(5));
    assert(upper(hull(Interval(2),Interval(5))) <= Boundary(6));

    // singleton
    assert(singleton(hull(Interval(2),Interval(2))) == true);
    assert(singleton(hull(Interval(2),Interval(3))) == false);

    // width
    assert(width(hull(Interval(2),Interval(2))) == Boundary(0));
    assert(width(hull(Interval(2),Interval(3))) == Boundary(1));   
}

template<class Interval_>
void generic_test_bigfloat_interval(){
    
    typedef typename CGAL::Bigfloat_interval_traits<Interval_>::Self BFIT;
    typedef typename BFIT::NT Interval;
    typedef typename BFIT::BF Boundary;
    
    const typename BFIT::Set_precision set_precsion = typename BFIT::Set_precision();
    const typename BFIT::Get_precision get_precsion = typename BFIT::Get_precision();
    const typename BFIT::Get_significant_bits get_significant_bits 
        = typename BFIT::Get_significant_bits();

    
    //TODO: move this into an new Interval_traits
    // get_significant_bits
    
    CGAL::set_precision(Interval(),15);
    assert(CGAL::get_precision(Interval())    == 15);
    assert(CGAL::set_precision(Interval(),23) == 15);
    assert(CGAL::set_precision(Interval(),70) == 23);

    //TODO: define what get_significant_bits should do and test is. Better name ?
    // just a compile check
    CGAL::get_significant_bits(Interval(3));
}



template<class Interval>
void generic_test_convert_to_bfi(){
    typedef typename CGAL::Get_arithmetic_kernel<Interval>::Arithmetic_kernel AK;
    typedef typename AK::Integer Integer; 
    typedef typename AK::Rational Rational; 
    typedef typename AK::Field_with_sqrt FWS;
    typedef CGAL::Sqrt_extension<Integer,Integer> EXT;
    
    assert(CGAL::convert_to_bfi(Integer(1)) == Interval(1));
    assert(CGAL::convert_to_bfi(Rational(1)) == Interval(1));
    assert(CGAL::convert_to_bfi(FWS(1)) == Interval(1));   
    assert(CGAL::convert_to_bfi(EXT(1)) == Interval(1));   
    
    typedef typename CGAL::Coercion_traits<Integer,Interval>::Type CT_Integer_type;
    BOOST_STATIC_ASSERT(( ::boost::is_same<CT_Integer_type, Interval>::value));
    typedef typename CGAL::Coercion_traits<Rational,Interval>::Type CT_Rational_type;
    BOOST_STATIC_ASSERT(( ::boost::is_same<CT_Rational_type, Interval>::value));
    typedef typename CGAL::Coercion_traits<FWS,Interval>::Type CT_FWS_type;
    BOOST_STATIC_ASSERT(( ::boost::is_same<CT_FWS_type, Interval>::value));
    typedef typename CGAL::Coercion_traits<EXT,Interval>::Type CT_EXT_type;
    BOOST_STATIC_ASSERT(( ::boost::is_same<CT_EXT_type, Interval>::value));
}

int main(){

#ifdef CGAL_USE_LEDA
    typedef CGAL::LEDA_arithmetic_kernel LEDA_AK;
    CGAL::test_interval<LEDA_AK::Bigfloat_interval>();
    generic_test_bigfloat_interval<LEDA_AK::Bigfloat_interval>();
    generic_test_convert_to_bfi<LEDA_AK::Bigfloat_interval>();
#endif 

#ifdef CGAL_USE_CORE
    typedef CGAL::CORE_arithmetic_kernel CORE_AK;
    CGAL::test_interval<CORE_AK::Bigfloat_interval>();
    generic_test_bigfloat_interval<CORE_AK::Bigfloat_interval>();
    generic_test_convert_to_bfi<CORE_AK::Bigfloat_interval>();
#endif 

}
