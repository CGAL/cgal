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

#ifndef CGAL_TEST_RATIONAL_TRAITS_H
#define CGAL_TEST_RATIONAL_TRAITS_H

CGAL_BEGIN_NAMESPACE

template <class Rational>
void test_rational_traits(){
    Rational x = Rational(7)/Rational(2);
    
    typedef Rational_traits<Rational> Rational_traits;
    typedef typename Rational_traits::RT RT;
    BOOST_STATIC_ASSERT((::boost::is_same<RT,RT>::value));
     
    CGAL_test_assert( Rational_traits().numerator(x) == RT(7));
    CGAL_test_assert( Rational_traits().denominator(x) == RT(2));
    CGAL_test_assert( Rational_traits().make_rational(RT(7),RT(2)) == x);
    CGAL_test_assert( Rational_traits().make_rational(x,x) == Rational(1));

    // gloabal function to_rational 
    x = CGAL::to_rational<Rational>(3.5);
    CGAL_test_assert( x == Rational(7)/Rational(2));
}

CGAL_END_NAMESPACE

#endif //  CGAL_TEST_RATIONAL_TRAITS_H
