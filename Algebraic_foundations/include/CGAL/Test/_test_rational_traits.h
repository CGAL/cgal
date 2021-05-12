// Copyright (c) 2006-2008 Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Michael Hemmer    <hemmer@mpi-inf.mpg.de>
//
// =============================================================================


// within this file FT ^= Fraction_traits<T>

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
