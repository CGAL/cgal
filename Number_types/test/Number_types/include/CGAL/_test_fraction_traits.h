// TODO: Add licence
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
// Author(s)     : Michael Hemmer  <mhemmer@uni-mainz.de>
//


// within this file FT ^= Fraction_traits<T>

#include <CGAL/basic.h>
#include <CGAL/Testsuite/assert.h>
#include <CGAL/to_rational.h>

#ifndef CGAL_TEST_FRACTION_TRAITS_H
#define CGAL_TEST_FRACTION_TRAITS_H

CGAL_BEGIN_NAMESPACE

template <class T>
void test_fraction_traits(){

    typedef CGAL::Fraction_traits<T> FT;    
    typedef typename FT::Type Type;
    typedef typename FT::Is_fraction  Is_fraction;
    typedef typename FT::Numerator_type  Num;
    typedef typename FT::Denominator_type Den;
    typedef typename FT::Common_factor Common_factor;
    typedef typename FT::Decompose Decompose;
    typedef typename FT::Compose   Compose;

    BOOST_STATIC_ASSERT( (::boost::is_same<Type,T>::value));
    BOOST_STATIC_ASSERT( (::boost::is_same<Is_fraction,Tag_true>::value));
    BOOST_STATIC_ASSERT(!(::boost::is_same<Common_factor,Null_functor>::value));
    BOOST_STATIC_ASSERT(!(::boost::is_same<Decompose,Null_functor>::value));
    BOOST_STATIC_ASSERT(!(::boost::is_same<Compose,Null_functor>::value));
    
    
    // Decompose
    Type frac = Type(7) / Type (5);
    Num num; 
    Den den;
    Decompose()(frac,num,den);
    CGAL_test_assert(num == Num(7));
    CGAL_test_assert(den == Num(5));
    CGAL_test_assert(frac == Compose()(num,den));
    
    // almost the same as gcd 
    Common_factor common_factor;
    CGAL_test_assert(common_factor(Den(0),Den(0)) == Den(0));
    CGAL_test_assert(common_factor(Den(1),Den(0)) == Den(1));
    CGAL_test_assert(common_factor(Den(-2),Den(0)) == Den(2));
    CGAL_test_assert(common_factor(Den(0),Den(-2)) == Den(2));
    CGAL_test_assert(common_factor(Den(12),Den(15)) == Den(3));
    CGAL_test_assert(common_factor(Den(-12),Den(15)) == Den(3));
    CGAL_test_assert(common_factor(Den(12),Den(-15)) == Den(3));
    CGAL_test_assert(common_factor(Den(-12),Den(-15)) == Den(3));    
}

CGAL_END_NAMESPACE

#endif //  CGAL_TEST_FRACTION_TRAITS_H
