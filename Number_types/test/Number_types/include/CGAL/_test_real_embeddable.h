// TODO: Add licence
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL: $
// $Id: $
// 
// Author(s)     : Lutz Kettner  <kettner@mpi-inf.mpg.de>
//                 Michael Hemmer <hemmer@mpi-inf.mpg.de>
//                 Michael Hemmer <mhemmer@uni-mainz.de>
//
// ============================================================================
//
//    \brief provides test functions for the \c RealEmbeddable concept of
//    number types.

#include <CGAL/basic.h>
#include <CGAL/_test_basic.h>

#include <cstddef>
#include <boost/static_assert.hpp>
#include <boost/type_traits.hpp>

#include <CGAL/Testsuite/assert.h>
#include <CGAL/tags.h>

//#include <CGAL/Algebraic_structure_traits.h>
#include <CGAL/Real_embeddable_traits.h>

#ifndef CGAL_TEST_REAL_EMBEDDABLE_H
#define CGAL_TEST_REAL_EMBEDDABLE_H

CGAL_BEGIN_NAMESPACE

    template<class Real_embeddable, class ToDouble>
    class Test_to_double {
    public:
        void operator() (ToDouble to_double) {
            typedef typename ToDouble::argument_type Argument_type;
            typedef typename ToDouble::result_type   Result_type;
            BOOST_STATIC_ASSERT(( ::boost::is_same<Real_embeddable, Argument_type>::value));  
            BOOST_STATIC_ASSERT(( ::boost::is_same<double, Result_type>::value));
            CGAL_test_assert(42.0 == to_double(Real_embeddable(42)));
        }
    };

    template<class Real_embeddable>
    class Test_to_double<Real_embeddable, CGAL::Null_tag> {
    public:
        void operator() (CGAL::Null_functor) {
            CGAL_error("To_double functor not implemented");
        }
    };
    
    template<class Real_embeddable, class To_interval>
    class Test_to_interval {
    public:
        void operator() (To_interval to_interval) {
            typedef typename To_interval::argument_type Argument_type;
            typedef typename To_interval::result_type   Result_type;
            typedef std::pair<double,double>  Interval_type;
            BOOST_STATIC_ASSERT(( ::boost::is_same<Real_embeddable, Argument_type>::value));
            BOOST_STATIC_ASSERT(( ::boost::is_same<Interval_type, Result_type>::value));

//            CGAL_test_assert(NiX::in(42.0,to_Interval(Real_embeddable(42))));
            // Instead of 'NiX::in':
            CGAL_test_assert( 42.0 >= to_interval( Real_embeddable(42) ).first );
            CGAL_test_assert( 42.0 <= to_interval( Real_embeddable(42) ).second );

            CGAL_test_assert(to_interval(Real_embeddable(42)).first > 41.99);
            CGAL_test_assert(to_interval(Real_embeddable(42)).second < 42.01);
            

	    /*
	    Real_embeddable notdouble = ipower(2,60);
            notdouble = notdouble + Real_embeddable(1);
            Interval test = to_Interval(notdouble);
            double lower = ipower(2.0,60);
            double upper = ipower(2.0,53);
            upper++;
            upper *= ipower(2.0,7);
	    std::cout << lower << "," << upper << std::endl;
	    std::cout << test.lower() << "," << test.upper() << std::endl;
            CGAL_test_assert( (in(lower,test) == true) && (in(upper,test) == true) );
            */
        }
    };
    
    template< class Real_embeddable >
    class Test_to_interval< Real_embeddable, CGAL::Null_tag> {
      public:
        void operator() (CGAL::Null_functor) {
            CGAL_assertion_msg(false, "To_Interval not implemented");
            // ok, nothing to test
        }
    };

template< class RE >
void test_ret_functor_arity() {
  typedef CGAL::Real_embeddable_traits<RE> RET;

  Test_functor_arity< typename RET::Is_finite >()(1);
  Test_functor_arity< typename RET::Abs >()(1);
  Test_functor_arity< typename RET::Compare >()(2);
  Test_functor_arity< typename RET::Sign >()(1);
  Test_functor_arity< typename RET::Is_zero >()(1);
  Test_functor_arity< typename RET::Is_positive >()(1);
  Test_functor_arity< typename RET::Is_negative >()(1);
  Test_functor_arity< typename RET::To_double >()(1);
  Test_functor_arity< typename RET::To_interval >()(1);
}


//! tests if \c Real_embeddable is a model for the \c RealComparable concept
//! and terminates the program with an error message if not.
template <class Real_embeddable>
void test_real_embeddable() {    
    typedef CGAL::Real_embeddable_traits<Real_embeddable> Traits;
    typedef typename Traits::Is_real_embeddable Is_real_embeddable;
    using CGAL::Tag_true;
    BOOST_STATIC_ASSERT(( ::boost::is_same< Is_real_embeddable, Tag_true>::value));

    test_ret_functor_arity< Real_embeddable >();
    typename Traits::Compare compare;
    typename Traits::Sign    sign;
    typename Traits::Abs     abs; 
    typename Traits::Is_finite  is_finite;
    typename Traits::Is_positive is_positive;
    typename Traits::Is_negative is_negative;
    typename Traits::Is_zero     is_zero;

    Real_embeddable a(-2);
    Real_embeddable b(1);
    Real_embeddable c(0);
    CGAL_test_assert( is_finite(a) );
    CGAL_test_assert( is_finite(b) );
    CGAL_test_assert( is_finite(c) );
    CGAL_test_assert( !is_positive(a) );
    CGAL_test_assert( is_positive(b) );
    CGAL_test_assert( !is_positive(c) );
    CGAL_test_assert( is_negative(a) );
    CGAL_test_assert( !is_negative(b) );
    CGAL_test_assert( !is_negative(c) );
    CGAL_test_assert( !is_zero(a) );
    CGAL_test_assert( !is_zero(b) );
    CGAL_test_assert( is_zero(c) );    
    CGAL_test_assert( a <  b);
    CGAL_test_assert( b >  a);
    CGAL_test_assert( a <= b);
    CGAL_test_assert( b >= a);
    CGAL_test_assert( a <= a);
    CGAL_test_assert( a >= a);
    CGAL_test_assert( compare(a,b) == CGAL::SMALLER);
    CGAL_test_assert( compare(b,a) == CGAL::LARGER);
    CGAL_test_assert( compare(a,a) == CGAL::EQUAL);
    CGAL_test_assert( sign(a) == CGAL::NEGATIVE);
    CGAL_test_assert( sign(b) == CGAL::POSITIVE);
    CGAL_test_assert( sign(c) == CGAL::ZERO);
    CGAL_test_assert( sign(a) <  sign(b));
    CGAL_test_assert( sign(b) >  sign(a));
    CGAL_test_assert( sign(a) <= sign(b));
    CGAL_test_assert( sign(b) >= sign(a));
    CGAL_test_assert( sign(a) <= sign(a));
    CGAL_test_assert( sign(a) >= sign(a));
    CGAL_test_assert( sign(c) <  sign(b));
    CGAL_test_assert( sign(b) >  sign(c));
    CGAL_test_assert( sign(c) <= sign(b));
    CGAL_test_assert( sign(b) >= sign(c)); 
    CGAL_test_assert( sign(c) <= sign(c));
    CGAL_test_assert( sign(c) >= sign(c));
    CGAL_test_assert( sign(a) <  sign(c));
    CGAL_test_assert( sign(c) >  sign(a));
    CGAL_test_assert( sign(a) <= sign(c));
    CGAL_test_assert( sign(c) >= sign(a));
    CGAL_test_assert( abs(a) == Real_embeddable(2));
    CGAL_test_assert( abs(b) == Real_embeddable(1));
    CGAL_test_assert( abs(c) == Real_embeddable(0));   
    
    // To_double --------------------------------------------------------------
    typename Traits::To_double  to_double;
    (void)to_double;
    Test_to_double<Real_embeddable, typename Traits::To_double> ttd;
    ttd(to_double);
    
    // To_Interval ------------------------------------------------------------
    typename Traits::To_interval  to_Interval;
    (void)to_Interval;
    Test_to_interval<Real_embeddable, typename Traits::To_interval> tti;
    tti(to_Interval);
    
    // additional functions     
    CGAL_test_assert( CGAL_NTS is_finite( Real_embeddable(1) ) );
    CGAL_test_assert( CGAL_NTS sign(Real_embeddable(-5))==CGAL::NEGATIVE);
    CGAL_test_assert( CGAL_NTS abs(Real_embeddable(-5))==Real_embeddable(5));
//    CGAL_test_assert(NiX::in(5.0,NiX::to_Interval(Real_embeddable(5))));
    CGAL_test_assert( CGAL_NTS compare(Real_embeddable(-5),Real_embeddable(6))==CGAL::SMALLER);
    CGAL_test_assert( CGAL_NTS is_positive(Real_embeddable(23)) );
    CGAL_test_assert( CGAL_NTS is_negative(Real_embeddable(-23)) );
    CGAL_test_assert( CGAL_NTS is_zero( Real_embeddable(0) ) );
    CGAL_test_assert( !CGAL_NTS is_zero( Real_embeddable(23) ) );

}

//! tests if \c Real_embeddable says it is not a model for the \c RealComparable 
//! concept and terminates the program with an error message if it 
//! actually is.
template <class Real_embeddable>
void test_not_real_embeddable() {
    typedef CGAL::Real_embeddable_traits<Real_embeddable> Traits;
    typedef typename Traits::Is_real_embeddable Is_real_embeddable;
    using CGAL::Tag_false;
    BOOST_STATIC_ASSERT(( ::boost::is_same< Is_real_embeddable, Tag_false>::value));
}


//template <class Real_embeddable, class CeilLog2Abs>
//void test_rounded_log2_abs(Real_embeddable zero, CGAL::Null_functor, CeilLog2Abs) {
//    typedef CGAL::Null_functor Null_functor;
//    BOOST_STATIC_ASSERT(( ::boost::is_same< CeilLog2Abs, Null_functor>::value));
//}
//
//template <class Real_embeddable, class FloorLog2Abs, class CeilLog2Abs>
//void test_rounded_log2_abs(Real_embeddable zero, FloorLog2Abs fl_log, CeilLog2Abs cl_log) {
//    typedef CGAL::Null_functor Null_functor;
//    BOOST_STATIC_ASSERT((!::boost::is_same< CeilLog2Abs, Null_functor>::value));
//
//    CGAL_test_assert( fl_log(Real_embeddable( 7)) == 2 );
//    CGAL_test_assert( cl_log(Real_embeddable( 7)) == 3 );
//    CGAL_test_assert( fl_log(Real_embeddable( 8)) == 3 );
//    CGAL_test_assert( cl_log(Real_embeddable( 8)) == 3 );
//    CGAL_test_assert( fl_log(Real_embeddable(-9)) == 3 );
//    CGAL_test_assert( cl_log(Real_embeddable(-9)) == 4 );
//
//    CGAL_test_assert( NiX::floor_log2_abs(Real_embeddable(  64)) == 6 );
//    CGAL_test_assert( NiX::ceil_log2_abs( Real_embeddable(  64)) == 6 );
//    CGAL_test_assert( NiX::floor_log2_abs(Real_embeddable(-126)) == 6 );
//    CGAL_test_assert( NiX::ceil_log2_abs( Real_embeddable(-126)) == 7 );
//}
//
//
////! tests that \c Floor_log2_abs and \c Ceil_log2_abs are
////! \c both CGAL::Null_functor or both working properly
////! (This is independent of the \c RealComparable concept)
//template <class Real_embeddable>
//void test_rounded_log2_abs() {
//
//    typedef typename NiX::Real_embeddable_traits<Real_embeddable>::Floor_log2_abs F;
//    typedef typename NiX::Real_embeddable_traits<Real_embeddable>::Ceil_log2_abs C;
//    test_rounded_log2_abs(Real_embeddable(0), F(), C());
//}

CGAL_END_NAMESPACE

#endif // CGAL_TEST_REAL_COMPARABLE_H
