// Copyright (c) 2006-2008 Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Lutz Kettner  <kettner@mpi-inf.mpg.de>
//                 Michael Hemmer <hemmer@mpi-inf.mpg.de>
//
// ============================================================================

//    \brief provides test functions for the \c RealEmbeddable concept of
//    number types.

#include <cstddef>
#include <CGAL/assertions.h>
#include <boost/type_traits.hpp>

#include <cassert>
#include <CGAL/tags.h>

//#include <CGAL/Algebraic_structure_traits.h>
#include <CGAL/Real_embeddable_traits.h>
#include <CGAL/Test/_test_algebraic_structure.h>

#define CGAL_SNAP_RET_FUNCTORS(RET)                      \
    typedef typename RET::Abs Abs;                       \
    typedef typename RET::Sgn Sgn;                       \
    typedef typename RET::Is_finite Is_finite;           \
    typedef typename RET::Is_positive Is_positive;       \
    typedef typename RET::Is_negative Is_negative;       \
    typedef typename RET::Is_zero Is_zero;               \
    typedef typename RET::Compare Compare;               \
    typedef typename RET::To_double To_double;           \
    typedef typename RET::To_interval To_interval;

#ifndef CGAL_TEST_REAL_EMBEDDABLE_H
#define CGAL_TEST_REAL_EMBEDDABLE_H

namespace CGAL {

    template<class Type, class ToDouble>
    class Test_to_double {
    public:
        void operator() (const ToDouble& to_double) {
            typedef typename ToDouble::argument_type Argument_type;
            typedef typename ToDouble::result_type   Result_type;
            CGAL_static_assertion(( ::boost::is_same<Type, Argument_type>::value));
            CGAL_USE_TYPE(Argument_type);
            CGAL_static_assertion(( ::boost::is_same<double, Result_type>::value));
            CGAL_USE_TYPE(Result_type);
            assert(42.0 == to_double(Type(42)));
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
            CGAL_static_assertion(( ::boost::is_same<Type, Argument_type>::value));
            CGAL_USE_TYPE(Argument_type);
            CGAL_static_assertion(( ::boost::is_same<Interval_type, Result_type>::value));
            CGAL_USE_TYPE(Result_type); CGAL_USE_TYPE(Interval_type);

//            assert(NiX::in(42.0,to_Interval(Type(42))));
            // Instead of 'NiX::in':
            assert( 42.0 >= to_interval( Type(42) ).first );
            assert( 42.0 <= to_interval( Type(42) ).second );

            assert(to_interval(Type(42)).first > 41.99);
            assert(to_interval(Type(42)).second < 42.01);

            // test neagtive numbers as well to catch obvious sign
            // errors
            assert( -42.0 >= to_interval( -Type(42) ).first );
            assert( -42.0 <= to_interval( -Type(42) ).second );

            assert(to_interval(-Type(42)).first < -41.99);
            assert(to_interval(-Type(42)).second > -42.01);


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
            assert( (in(lower,test) == true) && (in(upper,test) == true) );
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

// Several libraries do not want to enforce the use of std::min/max
// for the case that there is a special function for min/max
// However, this may be a problem for template CGAL::min/max in case NT is a
// CGAL type. These types have to overload CGAL::min/max.
template<typename NT>
void test_min_max(){
  using std:: min BOOST_PREVENT_MACRO_SUBSTITUTION ;
  using std:: max BOOST_PREVENT_MACRO_SUBSTITUTION ;
  NT x(1),y(2);
  assert(min BOOST_PREVENT_MACRO_SUBSTITUTION (x,y)==NT(1));
  assert(max BOOST_PREVENT_MACRO_SUBSTITUTION (x,y)==NT(2));
}

//! tests if \c Type is a model for the \c RealComparable concept
//! and terminates the program with an error message if not.
template <class Type>
void test_real_embeddable() {
  test_min_max<Type>();
    typedef CGAL::Real_embeddable_traits<Type> RET;
    CGAL_SNAP_RET_FUNCTORS(RET);
    typedef typename RET::Is_real_embeddable Is_real_embeddable;
    using CGAL::Tag_true;
    CGAL_static_assertion(( ::boost::is_same< Is_real_embeddable, Tag_true>::value));
    CGAL_USE_TYPE(Is_real_embeddable);

    typedef typename RET::Boolean Boolean;
    typedef typename RET::Sign Sign;
    typedef typename RET::Comparison_result Comparison_result;

    check_result_type(Is_finite()  ,Boolean());
    check_result_type(Is_positive(),Boolean());
    check_result_type(Is_negative(),Boolean());
    check_result_type(Is_zero(),Boolean());
    check_result_type(Sgn(),Sign());
    check_result_type(Compare(),Comparison_result());

    // test implicit interoperability with base type and proper default
    assert(Boolean()              == bool());
    assert(Sign()                 == CGAL::Sign());
    assert(Comparison_result()    == CGAL::Comparison_result());

    typename RET::Compare compare;
    const Sgn     sign = Sgn();
    const Abs     abs=Abs();
    const Is_finite  is_finite=Is_finite();
    const Is_positive is_positive=Is_positive();
    const Is_negative is_negative=Is_negative();
    const Is_zero     is_zero=Is_zero();

    Type a(-2);
    Type b(1);
    Type c(0);
    assert( is_finite(a) );
    assert( is_finite(b) );
    assert( is_finite(c) );
    assert( !is_positive(a) );
    assert( is_positive(b) );
    assert( !is_positive(c) );
    assert( is_negative(a) );
    assert( !is_negative(b) );
    assert( !is_negative(c) );
    assert( !is_zero(a) );
    assert( !is_zero(b) );
    assert( is_zero(c) );
    assert( a <  b);
    assert( b >  a);
    assert( a <= b);
    assert( b >= a);
    assert( a <= a);
    assert( a >= a);
    assert( compare(a,b) == CGAL::SMALLER);
    assert( compare(b,a) == CGAL::LARGER);
    assert( compare(a,a) == CGAL::EQUAL);
    assert( sign(a) == CGAL::NEGATIVE);
    assert( sign(b) == CGAL::POSITIVE);
    assert( sign(c) == CGAL::ZERO);
    assert( sign(a) <  sign(b));
    assert( sign(b) >  sign(a));
    assert( sign(a) <= sign(b));
    assert( sign(b) >= sign(a));
    assert( sign(a) <= sign(a));
    assert( sign(a) >= sign(a));
    assert( sign(c) <  sign(b));
    assert( sign(b) >  sign(c));
    assert( sign(c) <= sign(b));
    assert( sign(b) >= sign(c));
    assert( sign(c) <= sign(c));
    assert( sign(c) >= sign(c));
    assert( sign(a) <  sign(c));
    assert( sign(c) >  sign(a));
    assert( sign(a) <= sign(c));
    assert( sign(c) >= sign(a));
    assert( abs(a) == Type(2));
    assert( abs(b) == Type(1));
    assert( abs(c) == Type(0));

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
    assert( CGAL_NTS is_finite( Type(1) ) );
    assert( CGAL_NTS sign(Type(-5))==CGAL::NEGATIVE);
    assert( CGAL_NTS abs(Type(-5))==Type(5));
//    assert(NiX::in(5.0,NiX::to_interval(Type(5))));
    assert( CGAL_NTS compare(Type(-5),Type(6))==CGAL::SMALLER);
    assert( CGAL_NTS is_positive(Type(23)) );
    assert( CGAL_NTS is_negative(Type(-23)) );
    assert( CGAL_NTS is_zero( Type(0) ) );
    assert( !CGAL_NTS is_zero( Type(23) ) );

}

//! tests if \c Type says it is not a model for the \c RealComparable
//! concept and terminates the program with an error message if it
//! actually is.
template <class Type>
void test_not_real_embeddable() {
    typedef CGAL::Real_embeddable_traits<Type> RET;
    typedef typename RET::Is_real_embeddable Is_real_embeddable;
    using CGAL::Tag_false;
    CGAL_static_assertion(( ::boost::is_same< Is_real_embeddable, Tag_false>::value));
    CGAL_USE_TYPE(Is_real_embeddable);
}


//template <class Type, class CeilLog2Abs>
//void test_rounded_log2_abs(Type zero, CGAL::Null_functor, CeilLog2Abs) {
//    typedef CGAL::Null_functor Null_functor;
//    CGAL_static_assertion(( ::boost::is_same< CeilLog2Abs, Null_functor>::value));
//}
//
//template <class Type, class FloorLog2Abs, class CeilLog2Abs>
//void test_rounded_log2_abs(Type zero, FloorLog2Abs fl_log, CeilLog2Abs cl_log) {
//    typedef CGAL::Null_functor Null_functor;
//    CGAL_static_assertion((!::boost::is_same< CeilLog2Abs, Null_functor>::value));
//
//    assert( fl_log(Type( 7)) == 2 );
//    assert( cl_log(Type( 7)) == 3 );
//    assert( fl_log(Type( 8)) == 3 );
//    assert( cl_log(Type( 8)) == 3 );
//    assert( fl_log(Type(-9)) == 3 );
//    assert( cl_log(Type(-9)) == 4 );
//
//    assert( NiX::floor_log2_abs(Type(  64)) == 6 );
//    assert( NiX::ceil_log2_abs( Type(  64)) == 6 );
//    assert( NiX::floor_log2_abs(Type(-126)) == 6 );
//    assert( NiX::ceil_log2_abs( Type(-126)) == 7 );
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

} //namespace CGAL

#endif // CGAL_TEST_REAL_COMPARABLE_H
