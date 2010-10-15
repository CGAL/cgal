// ============================================================================
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
// $URL$
// $Id$
// author(s)     : Arno Eigenwillig <arno@mpi-inf.mpg.de>

#include <CGAL/Flattening_iterator.h>
#include <cassert>
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
        assert(*fi10 == i);
    }
    assert(fi10 == fi10_beyond);

    CGAL::Recursive_const_flattening<1, V2::const_iterator>
        ::Recursive_flattening_iterator fi21, fi21_beyond;
    fi21= CGAL::recursive_const_flattener<1>(const_cast<const V2&>(b1).end(), const_cast<const V2&>(b1).begin());
    fi21_beyond = CGAL::const_flattener(const_cast<const V2&>(b1).end(), const_cast<const V2&>(b1).end());
    for (i = 1; i <= 7; ++i, ++fi21) {
        assert(*fi21 == i);
    }
    assert(fi21 == fi21_beyond);

    CGAL::Recursive_const_flattening<2, V3::iterator>
        ::Recursive_flattening_iterator fi32, fi32_beyond;
    fi32 = CGAL::recursive_const_flattener<2>(c.end(),c.begin());
    fi32_beyond = CGAL::recursive_const_flattener<2>(c.end(),c.end());
    for (i = 1; i <= 17; ++i, ++fi32) {
        assert(*fi32 == i);
    }
    assert(fi32 == fi32_beyond);
}

int main() {
    test_const_flattening();
    return 0;
}

// EOF
