// Copyright (c) 2006-2008 Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
//
// Author(s)     : Lutz Kettner  <kettner@mpi-inf.mpg.de>
//                 Michael Hemmer <hemmer@mpi-inf.mpg.de>
//
// ============================================================================

//    \brief provides test functions for the \c RealEmbeddable concept of
//    number types.

#include <CGAL/basic.h>

#include <cstddef>
#include <CGAL/assertions.h>
#include <boost/type_traits.hpp>

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
            CGAL_assertion(42.0 == to_double(Type(42)));
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

//            CGAL_assertion(NiX::in(42.0,to_Interval(Type(42))));
            // Instead of 'NiX::in':
            CGAL_assertion( 42.0 >= to_interval( Type(42) ).first );
            CGAL_assertion( 42.0 <= to_interval( Type(42) ).second );

            CGAL_assertion(to_interval(Type(42)).first > 41.99);
            CGAL_assertion(to_interval(Type(42)).second < 42.01);
            
	    // test neagtive numbers as well to catch obvious sign
	    // errors
	    CGAL_assertion( -42.0 >= to_interval( -Type(42) ).first );
            CGAL_assertion( -42.0 <= to_interval( -Type(42) ).second );

            CGAL_assertion(to_interval(-Type(42)).first < -41.99);
            CGAL_assertion(to_interval(-Type(42)).second > -42.01);
     

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
            CGAL_assertion( (in(lower,test) == true) && (in(upper,test) == true) );
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
  CGAL_assertion(min BOOST_PREVENT_MACRO_SUBSTITUTION (x,y)==NT(1));
  CGAL_assertion(max BOOST_PREVENT_MACRO_SUBSTITUTION (x,y)==NT(2));
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
    CGAL_assertion(Boolean()              == bool());
    CGAL_assertion(Sign()                 == CGAL::Sign());
    CGAL_assertion(Comparison_result()    == CGAL::Comparison_result());
    
    CGAL_assertion_code(typename RET::Compare compare;)
    const Sgn     sign = Sgn();
    const Abs     abs=Abs(); 
    const Is_finite  is_finite=Is_finite();
    const Is_positive is_positive=Is_positive();
    const Is_negative is_negative=Is_negative();
    const Is_zero     is_zero=Is_zero();

    Type a(-2);
    Type b(1);
    Type c(0);
    CGAL_assertion( is_finite(a) );
    CGAL_assertion( is_finite(b) );
    CGAL_assertion( is_finite(c) );
    CGAL_assertion( !is_positive(a) );
    CGAL_assertion( is_positive(b) );
    CGAL_assertion( !is_positive(c) );
    CGAL_assertion( is_negative(a) );
    CGAL_assertion( !is_negative(b) );
    CGAL_assertion( !is_negative(c) );
    CGAL_assertion( !is_zero(a) );
    CGAL_assertion( !is_zero(b) );
    CGAL_assertion( is_zero(c) );    
    CGAL_assertion( a <  b);
    CGAL_assertion( b >  a);
    CGAL_assertion( a <= b);
    CGAL_assertion( b >= a);
    CGAL_assertion( a <= a);
    CGAL_assertion( a >= a);
    CGAL_assertion( compare(a,b) == CGAL::SMALLER);
    CGAL_assertion( compare(b,a) == CGAL::LARGER);
    CGAL_assertion( compare(a,a) == CGAL::EQUAL);
    CGAL_assertion( sign(a) == CGAL::NEGATIVE);
    CGAL_assertion( sign(b) == CGAL::POSITIVE);
    CGAL_assertion( sign(c) == CGAL::ZERO);
    CGAL_assertion( sign(a) <  sign(b));
    CGAL_assertion( sign(b) >  sign(a));
    CGAL_assertion( sign(a) <= sign(b));
    CGAL_assertion( sign(b) >= sign(a));
    CGAL_assertion( sign(a) <= sign(a));
    CGAL_assertion( sign(a) >= sign(a));
    CGAL_assertion( sign(c) <  sign(b));
    CGAL_assertion( sign(b) >  sign(c));
    CGAL_assertion( sign(c) <= sign(b));
    CGAL_assertion( sign(b) >= sign(c)); 
    CGAL_assertion( sign(c) <= sign(c));
    CGAL_assertion( sign(c) >= sign(c));
    CGAL_assertion( sign(a) <  sign(c));
    CGAL_assertion( sign(c) >  sign(a));
    CGAL_assertion( sign(a) <= sign(c));
    CGAL_assertion( sign(c) >= sign(a));
    CGAL_assertion( abs(a) == Type(2));
    CGAL_assertion( abs(b) == Type(1));
    CGAL_assertion( abs(c) == Type(0));   
    
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
    CGAL_assertion( CGAL_NTS is_finite( Type(1) ) );
    CGAL_assertion( CGAL_NTS sign(Type(-5))==CGAL::NEGATIVE);
    CGAL_assertion( CGAL_NTS abs(Type(-5))==Type(5));
//    CGAL_assertion(NiX::in(5.0,NiX::to_interval(Type(5))));
    CGAL_assertion( CGAL_NTS compare(Type(-5),Type(6))==CGAL::SMALLER);
    CGAL_assertion( CGAL_NTS is_positive(Type(23)) );
    CGAL_assertion( CGAL_NTS is_negative(Type(-23)) );
    CGAL_assertion( CGAL_NTS is_zero( Type(0) ) );
    CGAL_assertion( !CGAL_NTS is_zero( Type(23) ) );

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
//    CGAL_assertion( fl_log(Type( 7)) == 2 );
//    CGAL_assertion( cl_log(Type( 7)) == 3 );
//    CGAL_assertion( fl_log(Type( 8)) == 3 );
//    CGAL_assertion( cl_log(Type( 8)) == 3 );
//    CGAL_assertion( fl_log(Type(-9)) == 3 );
//    CGAL_assertion( cl_log(Type(-9)) == 4 );
//
//    CGAL_assertion( NiX::floor_log2_abs(Type(  64)) == 6 );
//    CGAL_assertion( NiX::ceil_log2_abs( Type(  64)) == 6 );
//    CGAL_assertion( NiX::floor_log2_abs(Type(-126)) == 6 );
//    CGAL_assertion( NiX::ceil_log2_abs( Type(-126)) == 7 );
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
