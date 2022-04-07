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
#include <CGAL/use.h>

#ifndef CGAL_TEST_FRACTION_TRAITS_H
#define CGAL_TEST_FRACTION_TRAITS_H

namespace CGAL {

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

    CGAL_USE_TYPE(Is_fraction);
    CGAL_static_assertion( (::boost::is_same<Type,T>::value));
    CGAL_static_assertion( (::boost::is_same<Is_fraction,Tag_true>::value));
    CGAL_static_assertion(!(::boost::is_same<Common_factor,Null_functor>::value));
    CGAL_static_assertion(!(::boost::is_same<Decompose,Null_functor>::value));
    CGAL_static_assertion(!(::boost::is_same<Compose,Null_functor>::value));


    // Decompose
    Type frac = Type(7) / Type (5);
    Num num;
    Den den;
    Decompose()(frac,num,den);
    assert(num == Num(7));
    assert(den == Num(5));
    assert(frac == Compose()(num,den));

    // almost the same as gcd
    Common_factor common_factor;
    assert(common_factor(Den(0),Den(0)) == Den(0));
    assert(common_factor(Den(1),Den(0)) == Den(1));
    assert(common_factor(Den(-2),Den(0)) == Den(2));
    assert(common_factor(Den(0),Den(-2)) == Den(2));
    assert(common_factor(Den(12),Den(15)) == Den(3));
    assert(common_factor(Den(-12),Den(15)) == Den(3));
    assert(common_factor(Den(12),Den(-15)) == Den(3));
    assert(common_factor(Den(-12),Den(-15)) == Den(3));
}

} //namespace CGAL

#endif //  CGAL_TEST_FRACTION_TRAITS_H
