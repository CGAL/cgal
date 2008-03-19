// ============================================================================
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
// $URL$
// $Id$
// Author(s)     : Lutz Kettner  <kettner@mpi-inf.mpg.de>
//                 Arno Eigenwillig <arno@mpi-inf.mpg.de>
//
// ============================================================================

#include <CGAL/basic.h>
#include <cassert>
#include <CGAL/algorithm.h>
#include <cstdlib>
#include <sstream>

int A[]  = {1,2,3,4,5};

int B1[] = {1,1,3};
int B2[] = {1,2,3};

CGAL::Comparison_result compare( int i, int j) {
    if ( i<j)
        return CGAL::SMALLER;
    if ( i>j)
        return CGAL::LARGER;
    return CGAL::EQUAL;
}

void test_lex_compare() {
    assert( CGAL::EQUAL   == CGAL::lexicographical_compare_three_valued( A, A+3, A, A+3, compare));
    assert( CGAL::SMALLER    == CGAL::lexicographical_compare_three_valued( A, A+3, A, A+4, compare));
    assert( CGAL::SMALLER    == CGAL::lexicographical_compare_three_valued( A, A+3, A, A+5, compare));
    assert( CGAL::LARGER == CGAL::lexicographical_compare_three_valued( A, A+4, A, A+3, compare));
    assert( CGAL::LARGER == CGAL::lexicographical_compare_three_valued( A, A+5, A, A+3, compare));

    assert( CGAL::SMALLER    == CGAL::lexicographical_compare_three_valued( B1, B1+3, B2, B2+3, compare));
    assert( CGAL::LARGER == CGAL::lexicographical_compare_three_valued( B2, B2+3, B1, B1+3, compare));
}

void test_output_range() {
    std::ostringstream os;
    std::ostream* sp;
    CGAL::set_ascii_mode(os);

    assert(os.str() == "");

    sp = &(CGAL::output_range(os, A, A, ":", "(", ")"));
    assert(os.str() == "");
    assert(sp == &os);
    os.str("");

    sp = &(CGAL::output_range(os, A, A+1, ":", "(", ")"));
    assert(os.str() == "(1)");
    assert(sp == &os);
    os.str("");

    sp = &(CGAL::output_range(os, A, A+3, ":", "(", ")"));
    assert(os.str() == "(1):(2):(3)");
    assert(sp == &os);
    os.str("");

    sp = &(CGAL::output_range(os, A, A+3));
    assert(os.str() == "1, 2, 3");
    assert(sp == &os);
}



int main() {
    test_lex_compare();
    test_output_range();

    return EXIT_SUCCESS;
}

// EOF
