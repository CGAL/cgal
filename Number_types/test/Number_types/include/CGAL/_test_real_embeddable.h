// TODO: Add licence
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
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

#define CGAL_SNAP_RET_FUNCTORS(RET)                      \
    typedef typename RET::Abs Abs;                       \
    typedef typename RET::Sign Sign;                     \
    typedef typename RET::Is_finite Is_finite;           \
    typedef typename RET::Is_positive Is_positive;       \
    typedef typename RET::Is_negative Is_negative;       \
    typedef typename RET::Is_zero Is_zero;               \
    typedef typename RET::Compare Compare;               \
    typedef typename RET::To_double To_double;           \
    typedef typename RET::To_interval To_interval;       

#ifndef CGAL_TEST_REAL_EMBEDDABLE_H
#define CGAL_TEST_REAL_EMBEDDABLE_H

CGAL_BEGIN_NAMESPACE

    template<class Type, class ToDouble>
    class Test_to_double {
    public:
        void operator() (const ToDouble& to_double) {
            typedef typename ToDouble::argument_type Argument_type;
            typedef typename ToDouble::result_type   Result_type;
            BOOST_STATIC_ASSERT(( ::boost::is_same<Type, Argument_type>::value));  
            BOOST_STATIC_ASSERT(( ::boost::is_same<double, Result_type>::value));
            CGAL_test_assert(42.0 == to_double(Type(42)));
        }
    };

    template<class Type>
    class Test_to_double<Type, CGAL::Null_tag> {
    public:
        void operator() (CGAL::Null_functor) {
            CGAL_error_msg("To_double functor not implemented");
        }
    };
    
    template<class Type, class To_interval>
    class Test_to_interval {
    public:
        void operator() (const To_interval& to_interval) {
            typedef typename To_interval::argument_type Argument_type;
            typedef typename To_interval::result_type   Result_type;
            typedef std::pair<double,double>  Interval_type;
            BOOST_STATIC_ASSERT(( ::boost::is_same<Type, Argument_type>::value));
            BOOST_STATIC_ASSERT(( ::boost::is_same<Interval_type, Result_type>::value));

//            CGAL_test_assert(NiX::in(42.0,to_Interval(Type(42))));
            // Instead of 'NiX::in':
            CGAL_test_assert( 42.0 >= to_interval( Type(42) ).first );
            CGAL_test_assert( 42.0 <= to_interval( Type(42) ).second );

            CGAL_test_assert(to_interval(Type(42)).first > 41.99);
            CGAL_test_assert(to_interval(Type(42)).second < 42.01);
            
	    // test neagtive numbers as well to catch obvious sign
	    // errors
	    CGAL_test_assert( -42.0 >= to_interval( -Type(42) ).first );
            CGAL_test_assert( -42.0 <= to_interval( -Type(42) ).second );

            CGAL_test_assert(to_interval(-Type(42)).first < -41.99);
            CGAL_test_assert(to_interval(-Type(42)).second > -42.01);
     

	    /*
	    Type notdouble = ipower(2,60);
            notdouble = notdouble + Type(1);
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
    
    template< class Type >
    class Test_to_interval< Type, CGAL::Null_tag> {
      public:
        void operator() (CGAL::Null_functor) {
            CGAL_error_msg( "To_Interval not implemented");
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


//! tests if \c Type is a model for the \c RealComparable concept
//! and terminates the program with an error message if not.
template <class Type>
void test_real_embeddable() {    
    typedef CGAL::Real_embeddable_traits<Type> RET;
    CGAL_SNAP_RET_FUNCTORS(RET);
    typedef typename RET::Is_real_embeddable Is_real_embeddable;
    using CGAL::Tag_true;
    BOOST_STATIC_ASSERT(( ::boost::is_same< Is_real_embeddable, Tag_true>::value));

    test_ret_functor_arity< Type >();
    typename RET::Compare compare;
    const Sign    sign = Sign();
    const Abs     abs=Abs(); 
    const Is_finite  is_finite=Is_finite();
    const Is_positive is_positive=Is_positive();
    const Is_negative is_negative=Is_negative();
    const Is_zero     is_zero=Is_zero();

    Type a(-2);
    Type b(1);
    Type c(0);
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
    CGAL_test_assert( abs(a) == Type(2));
    CGAL_test_assert( abs(b) == Type(1));
    CGAL_test_assert( abs(c) == Type(0));   
    
    // To_double --------------------------------------------------------------
    const To_double to_double = To_double();
    (void)to_double;
    Test_to_double<Type, To_double> ttd;
    ttd(to_double);
    
    // To_Interval ------------------------------------------------------------
    const To_interval to_interval = To_interval();
    (void)to_interval;
    Test_to_interval<Type, To_interval> tti;
    tti(to_interval);
    
    // additional functions     
    CGAL_test_assert( CGAL_NTS is_finite( Type(1) ) );
    CGAL_test_assert( CGAL_NTS sign(Type(-5))==CGAL::NEGATIVE);
    CGAL_test_assert( CGAL_NTS abs(Type(-5))==Type(5));
//    CGAL_test_assert(NiX::in(5.0,NiX::to_interval(Type(5))));
    CGAL_test_assert( CGAL_NTS compare(Type(-5),Type(6))==CGAL::SMALLER);
    CGAL_test_assert( CGAL_NTS is_positive(Type(23)) );
    CGAL_test_assert( CGAL_NTS is_negative(Type(-23)) );
    CGAL_test_assert( CGAL_NTS is_zero( Type(0) ) );
    CGAL_test_assert( !CGAL_NTS is_zero( Type(23) ) );

}

//! tests if \c Type says it is not a model for the \c RealComparable 
//! concept and terminates the program with an error message if it 
//! actually is.
template <class Type>
void test_not_real_embeddable() {
    typedef CGAL::Real_embeddable_traits<Type> RET;
    typedef typename RET::Is_real_embeddable Is_real_embeddable;
    using CGAL::Tag_false;
    BOOST_STATIC_ASSERT(( ::boost::is_same< Is_real_embeddable, Tag_false>::value));
}


//template <class Type, class CeilLog2Abs>
//void test_rounded_log2_abs(Type zero, CGAL::Null_functor, CeilLog2Abs) {
//    typedef CGAL::Null_functor Null_functor;
//    BOOST_STATIC_ASSERT(( ::boost::is_same< CeilLog2Abs, Null_functor>::value));
//}
//
//template <class Type, class FloorLog2Abs, class CeilLog2Abs>
//void test_rounded_log2_abs(Type zero, FloorLog2Abs fl_log, CeilLog2Abs cl_log) {
//    typedef CGAL::Null_functor Null_functor;
//    BOOST_STATIC_ASSERT((!::boost::is_same< CeilLog2Abs, Null_functor>::value));
//
//    CGAL_test_assert( fl_log(Type( 7)) == 2 );
//    CGAL_test_assert( cl_log(Type( 7)) == 3 );
//    CGAL_test_assert( fl_log(Type( 8)) == 3 );
//    CGAL_test_assert( cl_log(Type( 8)) == 3 );
//    CGAL_test_assert( fl_log(Type(-9)) == 3 );
//    CGAL_test_assert( cl_log(Type(-9)) == 4 );
//
//    CGAL_test_assert( NiX::floor_log2_abs(Type(  64)) == 6 );
//    CGAL_test_assert( NiX::ceil_log2_abs( Type(  64)) == 6 );
//    CGAL_test_assert( NiX::floor_log2_abs(Type(-126)) == 6 );
//    CGAL_test_assert( NiX::ceil_log2_abs( Type(-126)) == 7 );
//}
//
//
////! tests that \c Floor_log2_abs and \c Ceil_log2_abs are
////! \c both CGAL::Null_functor or both working properly
////! (This is independent of the \c RealComparable concept)
//template <class Type>
//void test_rounded_log2_abs() {
//
//    typedef typename NiX::Real_embeddable_traits<Type>::Floor_log2_abs F;
//    typedef typename NiX::Real_embeddable_traits<Type>::Ceil_log2_abs C;
//    test_rounded_log2_abs(Type(0), F(), C());
//}

CGAL_END_NAMESPACE

#endif // CGAL_TEST_REAL_COMPARABLE_H
