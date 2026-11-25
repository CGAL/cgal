// ============================================================================
//
// Copyright (c) 2001-2006 Max-Planck-Institut Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of EXACUS (http://www.mpi-inf.mpg.de/projects/EXACUS/).
//
// ----------------------------------------------------------------------------
//
// Library       : LiS
// File          : test/io.C
// LiS_release   : $Name:  $
// Revision      : $Revision$
// Revision_date : $Date$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Lutz Kettner  <kettner@mpi-inf.mpg.de>
//
// ============================================================================

#include <CGAL/IO/io.h>
#include <cassert>
#include <iostream>
#include <sstream>
#include <cstdlib>
#include <optional>


template <class NT>
void test_io(const NT& x){
    NT tmp;
    {
        std::ostringstream os;
        os << x;
        std::istringstream is(os.str());
        is >> tmp;
        assert( x == tmp );
    }{
        std::ostringstream os;
        os << ::CGAL::IO::oformat(x);
        std::istringstream is(os.str());
        is >> ::CGAL::IO::iformat(tmp);
        assert( x == tmp );
    }{
        std::ostringstream os;
        ::CGAL::write(os,x);
        std::istringstream is(os.str());
        ::CGAL::read(is,tmp);
        assert( x == tmp );
    }
}

// an empty optional is written as "--" and must be read back as such
template <class NT>
void test_io_optional(const NT& x){
    std::optional<NT> a, b(x), c(-x), d;
    std::ostringstream os;
    os << ::CGAL::IO::oformat(a) << ' ' << ::CGAL::IO::oformat(b) << ' '
       << ::CGAL::IO::oformat(c) << ' ' << ::CGAL::IO::oformat(d);
    std::istringstream is(os.str());
    std::optional<NT> ta(x), tb, tc, td(x);
    is >> ::CGAL::IO::iformat(ta) >> ::CGAL::IO::iformat(tb)
       >> ::CGAL::IO::iformat(tc) >> ::CGAL::IO::iformat(td);
    assert( is );
    assert( a == ta );
    assert( b == tb );
    assert( c == tc );
    assert( d == td );
}

int main() {
    static_assert(CGAL::Output_rep<int>::is_specialized == false);
#ifndef CGAL_DISABLE_FROM_CHARS
    static_assert(CGAL::Input_rep<int>::is_specialized == true);
 #else
    static_assert(CGAL::Input_rep<int>::is_specialized == false);
#endif
    std::cout << "test_io: short "<< std::endl;
    test_io<short>(12);
    test_io<short>(0);
    test_io<short>(-12);

    std::cout << "test_io: int "<< std::endl;
    test_io<int>(12);
    test_io<int>(0);
    test_io<int>(-12);

    std::cout << "test_io: long int "<< std::endl;
    test_io<long int>(12);
    test_io<long int>(0);
    test_io<long int>(-12);

    std::cout << "test_io: long long int "<< std::endl;
    test_io<long long int>(12);
    test_io<long long int>(0);
    test_io<long long int>(-12);

    std::cout << "test_io: float "<< std::endl;
    test_io<float>(12);
    test_io<float>(0);
    test_io<float>(-12);

    std::cout << "test_io: double "<< std::endl;
    test_io<double>(12);
    test_io<double>(0);
    test_io<double>(-12);

    std::cout << "test_io: long double "<< std::endl;
    test_io<long double>(12);
    test_io<long double>(0);
    test_io<long double>(-12);

    std::cout << "test_io: std::optional "<< std::endl;
    test_io_optional<int>(12);
    test_io_optional<double>(1.5);

    //TODO: add more
    std::cout << "ok" <<std::endl;

    return 0;
}
