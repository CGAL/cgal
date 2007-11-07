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

/*! \file include/NiX/test_real_comparable.h
    \brief provides test functions for the \c RealComparable concept of
    number types.
*/

#ifndef CGAL_TEST_REAL_COMPARABLE_H
#define CGAL_TEST_REAL_COMPARABLE_H

#include <CGAL/basic.h>
#include <CGAL/Testsuite/assert.h>
#include <CGAL/tags.h>
/*#include <NiX/basic.h>
#include <NiX/NT_traits.h>
#include <NiX/number_type_utils.h>*/
#include <cstddef>
#include <boost/static_assert.hpp>
#include <boost/type_traits.hpp>


CGAL_BEGIN_NAMESPACE

namespace CGALi {

    template<class NT, class ToDouble>
    class Test_to_double {
    public:
        void operator() (ToDouble to_double) {
            typedef typename ToDouble::argument_type Argument_type;
            typedef typename ToDouble::result_type   Result_type;
            BOOST_STATIC_ASSERT((::boost::is_same<NT, Argument_type>::value));
            BOOST_STATIC_ASSERT((::boost::is_same<double, Result_type>::value));
            CGAL_test_assert(42.0 == to_double(NT(42)));
        }
    };

    template<class NT>
    class Test_to_double<NT, ::CGAL::Null_functor> {
    public:
        void operator() (::CGAL::Null_functor) {
            CGAL_error_msg("To_double functor not implemented");
        }
    };
    
    template<class NT, class ToInterval>
    class Test_to_Interval {
    public:
        void operator() (ToInterval to_Interval) {
            typedef typename ToInterval::argument_type Argument_type;
            typedef typename ToInterval::result_type   Result_type;
            BOOST_STATIC_ASSERT((::boost::is_same<NT, Argument_type>::value));
            BOOST_STATIC_ASSERT((::boost::is_same< typename Argument_type::Interval, Result_type>::value));

            // TODO: NiX::in not available!?
            //CGAL_test_assert(NiX::in(42.0,to_Interval(NT(42))));
            CGAL_test_assert(to_Interval(NT(42)).lower() > 41.99);
            CGAL_test_assert(to_Interval(NT(42)).upper() < 42.01);

	    /*
	    NT notdouble = ipower(2,60);
            notdouble = notdouble + NT(1);
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
    
    template<class NT>
    class Test_to_Interval<NT, ::CGAL::Null_functor> {
    public:
        void operator() (::CGAL::Null_functor) {
            CGAL_error_msg("To_Interval not implemented");
            // ok, nothing to test
        }
    };


//! tests if \c NT is a model for the \c RealComparable concept
//! and terminates the program with an error message if not.
template <class NT>
void test_real_comparable() {    
    typedef CGAL::Real_embeddable_traits<NT> Traits;
    typedef typename Traits::Is_real_embeddable Is_real_comparable;
    using ::CGAL::Tag_true;
    BOOST_STATIC_ASSERT((::boost::is_same< Is_real_comparable, Tag_true>::value));
    typename Traits::Compare compare;
    typename Traits::Sign    sign;
    typename Traits::Abs     abs; 

    NT a(-2);
    NT b(1);
    NT c(0);
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
    CGAL_test_assert( abs(a) == NT(2));
    CGAL_test_assert( abs(b) == NT(1));
    CGAL_test_assert( abs(c) == NT(0));   
    
    // To_double --------------------------------------------------------------
    typename Traits::To_double  to_double;
    (void)to_double;
    typename CGALi::Test_to_double<NT, typename Traits::To_double> ttd;
    ttd(to_double);
    
    // To_Interval ------------------------------------------------------------
    typename Traits::To_Interval  to_Interval;
    (void)to_Interval;
    typename CGALi::Test_to_Interval<NT, typename Traits::To_Interval> tti;
    tti(to_Interval);
    
    // additional functions     
    CGAL_test_assert(CGAL::sign(NT(-5))==CGAL::NEGATIVE);
    CGAL_test_assert(CGAL::abs(NT(-5))==NT(5));
    // TODO: NiX::in not available!?
    //CGAL_test_assert(NiX::in(5.0,NiX::to_Interval(NT(5))));
    CGAL_test_assert(CGAL::compare(NT(-5),NT(6))==CGAL::SMALLER);

}

//! tests if \c NT says it is not a model for the \c RealComparable 
//! concept and terminates the program with an error message if it 
//! actually is.
template <class NT>
void test_not_real_comparable() {
    typedef CGAL::Real_embeddable_traits<NT> Traits;
    typedef typename Traits::Is_real_embeddable Is_real_comparable;
    using ::CGAL::Tag_false;
    BOOST_STATIC_ASSERT((::boost::is_same< Is_real_comparable, Tag_false>::value));
}


template <class NT, class CeilLog2Abs>
void test_rounded_log2_abs(NT zero, ::CGAL::Null_functor, CeilLog2Abs) {
    typedef ::CGAL::Null_functor Nulltype;
    BOOST_STATIC_ASSERT((::boost::is_same< CeilLog2Abs, Nulltype>::value));
}

template <class NT, class FloorLog2Abs, class CeilLog2Abs>
void test_rounded_log2_abs(NT zero, FloorLog2Abs fl_log, CeilLog2Abs cl_log) {
    typedef ::CGAL::Null_functor Null_functor;
    BOOST_STATIC_ASSERT((!::boost::is_same< CeilLog2Abs, Null_functor>::value));

    CGAL_test_assert( fl_log(NT( 7)) == 2 );
    CGAL_test_assert( cl_log(NT( 7)) == 3 );
    CGAL_test_assert( fl_log(NT( 8)) == 3 );
    CGAL_test_assert( cl_log(NT( 8)) == 3 );
    CGAL_test_assert( fl_log(NT(-9)) == 3 );
    CGAL_test_assert( cl_log(NT(-9)) == 4 );

    // TODO: floor_log2_abs etc. not available yet!?
    /*CGAL_test_assert( NiX::floor_log2_abs(NT(  64)) == 6 );
    CGAL_test_assert( NiX::ceil_log2_abs( NT(  64)) == 6 );
    CGAL_test_assert( NiX::floor_log2_abs(NT(-126)) == 6 );
    CGAL_test_assert( NiX::ceil_log2_abs( NT(-126)) == 7 );*/
}


//! tests that \c Floor_log2_abs and \c Ceil_log2_abs are
//! \c both ::CGAL::Null_functor or both working properly
//! (This is independent of the \c RealComparable concept)
template <class NT>
void test_rounded_log2_abs() {

    // TODO: floor_log2_abs etc. not available yet!?
    /*typedef typename NiX::NT_traits<NT>::Floor_log2_abs F;
    typedef typename NiX::NT_traits<NT>::Ceil_log2_abs C;
    test_rounded_log2_abs(NT(0), F(), C());*/
}

} // namespace CGALi


CGAL_END_NAMESPACE

#endif // CGAL_TEST_REAL_COMPARABLE_H
