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
// File          : test/algorithm_test.C
// LiS_release   : $Name:  $
// Revision      : $Revision$
// Revision_date : $Date$
//
// Author(s)     : Lutz Kettner  <kettner@mpi-inf.mpg.de>
//                 Arno Eigenwillig <arno@mpi-inf.mpg.de>
//
// ============================================================================

#include <CGAL/basic.h>
#include <CGAL/Testsuite/assert.h>
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
    CGAL_test_assert( CGAL::EQUAL   == CGAL::lexicographical_compare_three_valued( A, A+3, A, A+3, compare));
    CGAL_test_assert( CGAL::SMALLER    == CGAL::lexicographical_compare_three_valued( A, A+3, A, A+4, compare));
    CGAL_test_assert( CGAL::SMALLER    == CGAL::lexicographical_compare_three_valued( A, A+3, A, A+5, compare));
    CGAL_test_assert( CGAL::LARGER == CGAL::lexicographical_compare_three_valued( A, A+4, A, A+3, compare));
    CGAL_test_assert( CGAL::LARGER == CGAL::lexicographical_compare_three_valued( A, A+5, A, A+3, compare));

    CGAL_test_assert( CGAL::SMALLER    == CGAL::lexicographical_compare_three_valued( B1, B1+3, B2, B2+3, compare));
    CGAL_test_assert( CGAL::LARGER == CGAL::lexicographical_compare_three_valued( B2, B2+3, B1, B1+3, compare));
}

void test_output_range() {
    std::ostringstream os;
    std::ostream* sp;
    CGAL::set_ascii_mode(os);

    CGAL_test_assert(os.str() == "");

    sp = &(CGAL::output_range(os, A, A, ":", "(", ")"));
    CGAL_test_assert(os.str() == "");
    CGAL_test_assert(sp == &os);
    os.str("");
    
    sp = &(CGAL::output_range(os, A, A+1, ":", "(", ")"));
    CGAL_test_assert(os.str() == "(1)");
    CGAL_test_assert(sp == &os);
    os.str("");

    sp = &(CGAL::output_range(os, A, A+3, ":", "(", ")"));
    CGAL_test_assert(os.str() == "(1):(2):(3)");
    CGAL_test_assert(sp == &os);
    os.str("");

    sp = &(CGAL::output_range(os, A, A+3));
    CGAL_test_assert(os.str() == "1, 2, 3");
    CGAL_test_assert(sp == &os);
}



int main() {
    test_lex_compare();
    test_output_range();
    
    return EXIT_SUCCESS;
}

// EOF
