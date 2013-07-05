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
//
// Author(s)     : Michael Hemmer    <hemmer@mpi-inf.mpg.de>
//
// =============================================================================


// within this file FT ^= Fraction_traits<T>

#include <CGAL/basic.h>
#include <cassert>
#include <CGAL/to_rational.h>

#ifndef CGAL_TEST_RATIONAL_TRAITS_H
#define CGAL_TEST_RATIONAL_TRAITS_H

namespace CGAL {

template <class Rational>
void test_rational_traits(){
    Rational x = Rational(7)/Rational(2);
    
    typedef Rational_traits<Rational> Rational_traits;
    typedef typename Rational_traits::RT RT;
    CGAL_static_assertion((::boost::is_same<RT,RT>::value));
     
    assert( Rational_traits().numerator(x) == RT(7));
    assert( Rational_traits().denominator(x) == RT(2));
    assert( Rational_traits().make_rational(RT(7),RT(2)) == x);
    assert( Rational_traits().make_rational(x,x) == Rational(1));
    assert( Rational_traits().make_rational(x) == x);
    assert( Rational_traits().make_rational(std::make_pair(x,x)) == Rational(1));
    assert( Rational_traits().make_rational(std::make_pair(7,RT(2))) == x);

    // gloabal function to_rational 
    x = CGAL::to_rational<Rational>(3.5);
    assert( x == Rational(7)/Rational(2));
}

} //namespace CGAL

#endif //  CGAL_TEST_RATIONAL_TRAITS_H
