// Copyright (c) 2004  Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Lutz Kettner  <kettner@mpi-sb.mpg.de>
//                 Andreas Meyer <ameyer@mpi-sb.mpg.de>

// enable invariant checking
#define SEGMENT_TREE_CHECK_INVARIANTS 1

#include <CGAL/box_intersection_d.h>
#include <CGAL/Timer.h>
#include <iostream>
#include <sstream>
#include <iterator>
#include <cstdio>

#include "util.h"

bool test_failed = false;

template< class NT, unsigned int DIM, bool CLOSED >
struct _test {
    typedef Util< NT, DIM, CLOSED > Util;


void
operator()( const char* filename1, const char* filename2 )
{
    typename Util::Box_container boxes1, boxes2;
    typename Util::Result_container result_all_pairs, result_tree;
    FILE *infile1, *infile2;
    infile1 = fopen( filename1, "r");
    infile2 = fopen( filename2, "r");

    Util::readBoxesFromFile( infile1, boxes1 );
    Util::readBoxesFromFile( infile2, boxes2 );

    std::cout << std::endl;
    typename Util::Storage_callback<>
        callback1( result_all_pairs ),
        callback2( result_tree );

    std::cout << "all pairs ...... " << std::flush;
    CGAL::Timer timer;
    timer.start();
    CGAL::Box_intersection_d::all_pairs( boxes1.begin(), boxes1.end(),
                                         boxes2.begin(), boxes2.end(),
                                         callback1, Util::Traits(), 2 );
    timer.stop();
    std::cout << "got " << callback1.counter << " intersections in "
              << timer.time() << " seconds." << std::endl;


    std::cout << "one way scan ...... " << std::flush;
    timer.reset();
    timer.start();
    CGAL::Box_intersection_d::one_way_scan( boxes1.begin(), boxes1.end(),
                                            boxes2.begin(), boxes2.end(),
                                            callback2, Util::Traits(), 2 );
    CGAL::Box_intersection_d::one_way_scan( boxes2.begin(), boxes2.end(),
                                            boxes1.begin(), boxes1.end(),
                                            callback2, Util::Traits(), 2 );
    timer.stop();
    std::cout << "got " << callback2.counter << " intersections in "
              << timer.time() << " seconds." << std::endl;
    callback2.counter = 0;
    result_tree.clear();

    std::cout << "segment tree ... " << std::flush;
    timer.reset();
    timer.start();
    const unsigned int n = boxes1.size();
    const unsigned int cutoff = n < 2000 ? 6 : n / 100;
    CGAL::box_intersection_custom_predicates_d( boxes1.begin(), boxes1.end(),
                                                boxes2.begin(), boxes2.end(),
                                                callback2, Util::Traits(), cutoff );
    timer.stop();
    std::cout << "got " << callback2.counter << " intersections in "
              << timer.time() << " seconds." << std::endl;

    if( callback1.counter != callback2.counter ) {
        unsigned int missing    = Util::countMissingItems( result_all_pairs,
                                                     result_tree );
        unsigned int duplicates = Util::countDuplicates( result_tree );
        std::cout << "!! failed !! " << missing  << " missing and "
             << duplicates << " duplicate intersections in tree result."
             << std::endl;
        test_failed = true;
    }
    else
        std::cout << "--- passed --- " << std::endl;
    fclose( infile1 );
    fclose( infile2 );
}

}; // end class test

int main( int argc, char ** argv ) {
    _test<int,3,true> test;
    for( unsigned int n = 1; n <= 6; ++n ) {
        std::stringstream file1, file2;
        file1 << "data/test" << n << "_set1.box" << std::ends;
        file2 << "data/test" << n << "_set2.box" << std::ends;
        test( file1.str().c_str(), file2.str().c_str() );
    }
    if ( test_failed)
        return 1;
    return 0;
}

