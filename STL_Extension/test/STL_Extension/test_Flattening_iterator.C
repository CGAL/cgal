// ============================================================================
//
// Copyright (c) 2001-2006 Max-Planck-Institut Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of EXACUS (http://www.mpi-inf.mpg.de/projects/EXACUS/);
// you may redistribute it under the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with EXACUS.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// ----------------------------------------------------------------------------
//
// Library       : LiS
// File          : test/Flattening_iterator.C
// LiS_release   : $Name:  $
// Revision      : $Revision$
// Revision_date : $Date$
//
// Author(s)     : Arno Eigenwillig <arno@mpi-inf.mpg.de>
//
// ============================================================================

#include <CGAL/Flattening_iterator.h>
#include <CGAL/Testsuite/assert.h>
#include <vector>
#include <functional>

void test_const_flattening() {
    typedef std::vector<int> V1;
    typedef std::vector<V1>  V2;
    typedef std::vector<V2>  V3;

    V1 a1, a2, a3, a4, a5, a6;

    a1.push_back(1); a1.push_back(2); a1.push_back(3); a1.push_back(4);
    a2.push_back(5); a2.push_back(6); a2.push_back(7);
    a3.push_back(8);
    a4.push_back(9); a4.push_back(10); a4.push_back(11); a4.push_back(12);
    a5.push_back(13); a5.push_back(14);
    a6.push_back(15); a6.push_back(16); a6.push_back(17);

    V2 b1, b2, b3;

    b1.push_back(a1); b1.push_back(a2);
    b2.push_back(a3);
    b3.push_back(a4); b3.push_back(a5); b3.push_back(a6);

    V3 c;
    c.push_back(b1); c.push_back(b2); c.push_back(b3);

    int i;

    CGAL::Recursive_const_flattening<0, V1::const_iterator>
        ::Recursive_flattening_iterator fi10, fi10_beyond;
    fi10 = CGAL::recursive_const_flattener<0>( a1.end(),a1.begin() );
    fi10_beyond = CGAL::recursive_const_flattener<0>( a1.end(),a1.end() );
    for (i = 1; i <= 4; ++i, ++fi10) {
        CGAL_test_assert(*fi10 == i);
    }
    CGAL_test_assert(fi10 == fi10_beyond);

    CGAL::Recursive_const_flattening<1, V2::const_iterator>
        ::Recursive_flattening_iterator fi21, fi21_beyond;
    fi21= CGAL::recursive_const_flattener<1>(const_cast<const V2&>(b1).end(), const_cast<const V2&>(b1).begin());
    fi21_beyond = CGAL::const_flattener(const_cast<const V2&>(b1).end(), const_cast<const V2&>(b1).end());
    for (i = 1; i <= 7; ++i, ++fi21) {
        CGAL_test_assert(*fi21 == i);
    }
    CGAL_test_assert(fi21 == fi21_beyond);

    CGAL::Recursive_const_flattening<2, V3::iterator>
        ::Recursive_flattening_iterator fi32, fi32_beyond;
    fi32 = CGAL::recursive_const_flattener<2>(c.end(),c.begin());
    fi32_beyond = CGAL::recursive_const_flattener<2>(c.end(),c.end());
    for (i = 1; i <= 17; ++i, ++fi32) {
        CGAL_test_assert(*fi32 == i);
    }
    CGAL_test_assert(fi32 == fi32_beyond);
}

int main(int argc, char** argv) {
    test_const_flattening();
    return 0;
}

// EOF
