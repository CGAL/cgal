// Copyright (c) 2006-2007 Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Michael Hemmer, Dominik Hülse
//
// ============================================================================

#include <CGAL/Arithmetic_kernel.h>
#include <CGAL/Algebraic_structure_traits.h>
#include <cassert>
#include <CGAL/extended_euclidean_algorithm.h>

//#define WITH_OUTPUT 1

template<class Integer>
void test_extended_euclidean_algorithm() {
    Integer a, b, e, u, v;

    // common factor is 2
    a = Integer(1008);
    b = Integer(4352);
    e = CGAL::extended_euclidean_algorithm(a, b,u ,v );
    assert(Integer(16) == e);
    assert(e == a*u + b*v);


    a = Integer(-1008);
    b = Integer(4352);
    e = CGAL::extended_euclidean_algorithm(a, b,u ,v );
    assert(Integer(16) == e);
    assert(e == a*u + b*v);

    a = Integer(1008);
    b = Integer(-4352);
    e = CGAL::extended_euclidean_algorithm(a, b,u ,v );
    assert(Integer(16) == e);
    assert(e == a*u + b*v);

    a = Integer(-1008);
    b = Integer(-4352);
    e = CGAL::extended_euclidean_algorithm(a, b,u ,v );
    assert(Integer(16) == e);
    assert(e == a*u + b*v);

    a = CGAL::integral_division(a,e);
    b = CGAL::integral_division(b,e);
    e = CGAL::extended_euclidean_algorithm(a, b,u ,v );
    assert(Integer(1) == e);
    assert(e == a*u + b*v);

    // special cases
    // common factor is 1
    a = Integer(17);
    b = Integer(13);
    e = CGAL::extended_euclidean_algorithm(a, b,u ,v );
    assert(Integer(1) == e);
    assert(e == a*u + b*v);

    // one number is 0
    a = Integer(0);
    b = Integer(13);
    e = CGAL::extended_euclidean_algorithm(a, b,u ,v );
    assert(b == e);
    assert(e == a*u + b*v);

    a = Integer(24);
    b = Integer(0);
    e = CGAL::extended_euclidean_algorithm(a, b,u ,v );
    assert(a == e);
    assert(e == a*u + b*v);

    // both numbers are 0
    a = Integer(0);
    b = Integer(0);
    e = CGAL::extended_euclidean_algorithm(a, b,u ,v );
    assert(b == e);
    assert(e == a*u + b*v);
}

int main(){
    test_extended_euclidean_algorithm<long>();
    return 0;
}


// EOF
