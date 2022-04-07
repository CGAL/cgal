// Copyright (c) 2004  Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Lutz Kettner  <kettner@mpi-sb.mpg.de>
//                 Andreas Meyer <ameyer@mpi-sb.mpg.de>

// enable invariant checking
#define CGAL_SEGMENT_TREE_CHECK_INVARIANTS 1

#ifndef _CRT_SECURE_NO_DEPRECATE
#define _CRT_SECURE_NO_DEPRECATE 1
#endif

#ifndef _CRT_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_WARNINGS 1
#endif

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
    typedef Util< NT, DIM, CLOSED > Uti1;


void
operator()( const char* filename1, const char* filename2 )
{
    typename Uti1::Box_container boxes1, boxes2;
    typename Uti1::Result_container result_all_pairs, result_scan, result_tree;
    std::FILE *infile1, *infile2;
    infile1 = std::fopen( filename1, "r");
    infile2 = std::fopen( filename2, "r");

    Uti1::readBoxesFromFile( infile1, boxes1 );
    Uti1::readBoxesFromFile( infile2, boxes2 );

    std::cout << std::endl;
    unsigned int c0 = 0, c1 = 0, c2 = 0, c3 = 0;
    typename Uti1::Counter_callback callback0(c0);
    typename Uti1::template Storage_callback<typename Uti1::Result_container>
      callback1( result_all_pairs, c1 ),
      callback2( result_scan, c2 ),
      callback3( result_tree, c3 );

    // invoke each interface routine at least once, to check if it still
    // compiles
    CGAL::box_intersection_all_pairs_d(
                boxes1.begin(), boxes1.end(),
                boxes2.begin(), boxes2.end(),
                callback0,
                CLOSED ?
                   CGAL::Box_intersection_d::CLOSED :
                   CGAL::Box_intersection_d::HALF_OPEN );
    std::cout << "all pairs ......... " << std::flush;
    CGAL::Timer timer;
    timer.start();
    CGAL::Box_intersection_d::all_pairs(
                boxes1.begin(), boxes1.end(),
                boxes2.begin(), boxes2.end(),
                callback1, typename Uti1::Traits() );
    timer.stop();
    std::cout << "got " << callback1.get_counter() << " intersections in "
              << timer.time() << " seconds." << std::endl;


    std::cout << "one way scan ...... " << std::flush;
    timer.reset();
    timer.start();
    CGAL::Box_intersection_d::one_way_scan( boxes1.begin(), boxes1.end(),
                                            boxes2.begin(), boxes2.end(),
                                            callback2,
                                            typename Uti1::Traits(), 2 );
    CGAL::Box_intersection_d::one_way_scan( boxes2.begin(), boxes2.end(),
                                            boxes1.begin(), boxes1.end(),
                                            callback2,
                                            typename Uti1::Traits(), 2 );
    timer.stop();
    std::cout << "got " << callback2.get_counter() << " intersections in "
              << timer.time() << " seconds." << std::endl;

    std::cout << "segment tree ...... " << std::flush;
    timer.reset();
    timer.start();
    const std::size_t n = boxes1.size();
    const std::size_t cutoff = n < 2000 ? 6 : n / 100;
    CGAL::box_intersection_custom_predicates_d( boxes1.begin(), boxes1.end(),
                                                boxes2.begin(), boxes2.end(),
                                                callback3,
                                                typename Uti1::Traits(),
                                                cutoff );
    timer.stop();
    std::cout << "got " << callback3.get_counter() << " intersections in "
              << timer.time() << " seconds." << std::endl;

    if( callback1.get_counter() != callback2.get_counter() ||
        callback1.get_counter() != callback3.get_counter() )
    {
        unsigned int missing    = Uti1::countMissingItems( result_all_pairs,
                                                           result_tree );
        unsigned int duplicates = Uti1::countDuplicates( result_tree );
        std::cout << "!! failed !! " << missing  << " missing and "
             << duplicates << " duplicate intersections in tree result."
             << std::endl;
        test_failed = true;
    }
    else
        std::cout << "--- passed --- " << std::endl;
    std::fclose( infile1 );
    std::fclose( infile2 );
}

}; // end class test

int main() {
    _test<int,3,true> test1;
    _test<int,3,false> test2;
     for( unsigned int n = 1; n <= 6; ++n ) {
        std::stringstream file1, file2;
        file1 << "data/test" << n << "_set1.box" << std::ends;
        file2 << "data/test" << n << "_set2.box" << std::ends;
        test1( file1.str().c_str(), file2.str().c_str() );
        test2( file1.str().c_str(), file2.str().c_str() );
    }
    if ( test_failed)
        return 1;
    return 0;
}

