// Copyright (c) 2006-2009 Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     :  Michael Hemmer <hemmer@mpi-inf.mpg.de>
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

#include <cassert>
#include <CGAL/tags.h>
#include <cstddef>
#include <CGAL/assertions.h>
#include <boost/type_traits.hpp>


namespace CGAL {

namespace internal {

    template<class NT, class ToDouble>
    class Test_to_double {
    public:
        void operator() (ToDouble to_double) {
            typedef typename ToDouble::argument_type Argument_type;
            typedef typename ToDouble::result_type   Result_type;
            CGAL_static_assertion((::boost::is_same<NT, Argument_type>::value));
            CGAL_static_assertion((::boost::is_same<double, Result_type>::value));
            assert(42.0 == to_double(NT(42)));
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
            CGAL_static_assertion((::boost::is_same<NT, Argument_type>::value));
            CGAL_static_assertion((::boost::is_same< typename Argument_type::Interval, Result_type>::value));

            // TODO: NiX::in not available!?
            //assert(NiX::in(42.0,to_Interval(NT(42))));
            assert(to_Interval(NT(42)).lower() > 41.99);
            assert(to_Interval(NT(42)).upper() < 42.01);

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
            assert( (in(lower,test) == true) && (in(upper,test) == true) );
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
    CGAL_static_assertion((::boost::is_same< Is_real_comparable, Tag_true>::value));
    typename Traits::Compare compare;
    typename Traits::Sign    sign;
    typename Traits::Abs     abs;

    NT a(-2);
    NT b(1);
    NT c(0);
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
    assert( abs(a) == NT(2));
    assert( abs(b) == NT(1));
    assert( abs(c) == NT(0));

    // To_double --------------------------------------------------------------
    typename Traits::To_double  to_double;
    (void)to_double;
    typename internal::Test_to_double<NT, typename Traits::To_double> ttd;
    ttd(to_double);

    // To_Interval ------------------------------------------------------------
    typename Traits::To_Interval  to_Interval;
    (void)to_Interval;
    typename internal::Test_to_Interval<NT, typename Traits::To_Interval> tti;
    tti(to_Interval);

    // additional functions
    assert(CGAL::sign(NT(-5))==CGAL::NEGATIVE);
    assert(CGAL::abs(NT(-5))==NT(5));
    // TODO: NiX::in not available!?
    //assert(NiX::in(5.0,NiX::to_Interval(NT(5))));
    assert(CGAL::compare(NT(-5),NT(6))==CGAL::SMALLER);

}

//! tests if \c NT says it is not a model for the \c RealComparable
//! concept and terminates the program with an error message if it
//! actually is.
template <class NT>
void test_not_real_comparable() {
    typedef CGAL::Real_embeddable_traits<NT> Traits;
    typedef typename Traits::Is_real_embeddable Is_real_comparable;
    using ::CGAL::Tag_false;
    CGAL_static_assertion((::boost::is_same< Is_real_comparable, Tag_false>::value));
}


template <class NT, class CeilLog2Abs>
void test_rounded_log2_abs(NT zero, ::CGAL::Null_functor, CeilLog2Abs) {
    typedef ::CGAL::Null_functor Nulltype;
    CGAL_static_assertion((::boost::is_same< CeilLog2Abs, Nulltype>::value));
}

template <class NT, class FloorLog2Abs, class CeilLog2Abs>
void test_rounded_log2_abs(NT zero, FloorLog2Abs fl_log, CeilLog2Abs cl_log) {
    typedef ::CGAL::Null_functor Null_functor;
    CGAL_static_assertion((!::boost::is_same< CeilLog2Abs, Null_functor>::value));

    assert( fl_log(NT( 7)) == 2 );
    assert( cl_log(NT( 7)) == 3 );
    assert( fl_log(NT( 8)) == 3 );
    assert( cl_log(NT( 8)) == 3 );
    assert( fl_log(NT(-9)) == 3 );
    assert( cl_log(NT(-9)) == 4 );

    // TODO: floor_log2_abs etc. not available yet!?
    /*assert( NiX::floor_log2_abs(NT(  64)) == 6 );
    assert( NiX::ceil_log2_abs( NT(  64)) == 6 );
    assert( NiX::floor_log2_abs(NT(-126)) == 6 );
    assert( NiX::ceil_log2_abs( NT(-126)) == 7 );*/
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

} // namespace internal


} //namespace CGAL

#endif // CGAL_TEST_REAL_COMPARABLE_H
