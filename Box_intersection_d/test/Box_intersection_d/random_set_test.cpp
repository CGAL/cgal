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

#include <CGAL/box_intersection_d.h>
#include <CGAL/Timer.h>
#include <iostream>
#include <iterator>

#include "util.h"

static unsigned int failed = 0;

template< class NT, std::size_t DIM, bool CLOSED >
struct _test {
typedef Util< NT, DIM, CLOSED > Uti1;

static void
test_n( unsigned int n,
        CGAL::Box_intersection_d::Setting
                             setting = CGAL::Box_intersection_d::BIPARTITE )
{
    typename Uti1::Box_container boxes1, boxes2;
    const unsigned int allpairs_max = 5000;
    //Result_container result_allpairs, result_scanner, result_tree;
    std::cout << "generating random box sets with size " << n
              << " ... " << std::flush;
    Uti1::fill_boxes( n, boxes1 );
    bool bipartite = setting == CGAL::Box_intersection_d::BIPARTITE;

    if( bipartite )
        Uti1::fill_boxes( n, boxes2 );
    else
        boxes2 = boxes1;
    std::cout << std::endl;
    unsigned int c0 = 0, c1 = 0, c2 = 0;
    typename Uti1::Counter_callback callback0(c0), callback1(c1), callback2(c2);
    CGAL::Timer timer;

    if( n < allpairs_max ) {
        std::cout << "all pairs ... " << std::flush;
        timer.start();
        CGAL::Box_intersection_d::all_pairs(
                        boxes1.begin(), boxes1.end(),
                        boxes2.begin(), boxes2.end(),
                        callback0, typename Uti1::Traits(), bipartite == false );
        timer.stop();
        std::cout << "got " << callback0.get_counter() << " intersections in "
                  << timer.time() << " seconds."
            << std::endl;
        timer.reset();
    }
    std::cout << "one way scan ... " << std::flush;
    timer.start();
    CGAL::Box_intersection_d::one_way_scan( boxes1.begin(), boxes1.end(),
                                            boxes2.begin(), boxes2.end(),
                                            callback1, typename Uti1::Traits(),
                                            DIM - 1 );
    if( bipartite )
        CGAL::Box_intersection_d::one_way_scan( boxes2.begin(), boxes2.end(),
                                                boxes1.begin(), boxes1.end(),
                                                callback1,
                                                typename Uti1::Traits(),
                                                DIM - 1);
    timer.stop();
    std::cout << "got " << callback1.get_counter() << " intersections in "
              << timer.time() << " seconds."
              << std::endl;

    std::cout << "segment tree ... " << std::flush;
    timer.reset();
    timer.start();
    unsigned int cutoff = n < 200 ? 6 : n < 2000 ? 20 : n / 50;
    CGAL::box_intersection_custom_predicates_d( boxes1.begin(), boxes1.end(),
                                      boxes2.begin(), boxes2.end(),
                                      callback2, typename Uti1::Traits(),
                                      cutoff, setting );
    timer.stop();
    std::cout << "got " << callback2.get_counter() << " intersections in "
              << timer.time() << " seconds." << std::endl;

    if( callback1.get_counter() != callback2.get_counter() ||
        ( n < allpairs_max && callback0.get_counter() != callback1.get_counter() ) )
    {
        ++failed;
        std::cout << "!! failed !! " << std::endl;
    } else
        std::cout << "--- passed --- " << std::endl;
}

void operator()() {
    std::cout << "-------------------------" << std::endl;
    std::cout << "DIM = " << DIM << std::endl;
    std::cout << "-------------------------" << std::endl;
    for( unsigned int n = 4; n < 20000; n *= 6) {
        std::cout << "bipartite case: " << std::endl;
        test_n( n, CGAL::Box_intersection_d::BIPARTITE );
        std::cout << "complete case: " << std::endl;
        test_n( n, CGAL::Box_intersection_d::COMPLETE );
    }
}

}; // end struct  test

template< class NT >
struct test
{
    _test< NT, 1, true >   t0;
    //_test< NT, 2, true >   t1;
    _test< NT, 3, true >   t2;
    _test< NT, 4, true >   t3;
    //_test< NT, 10, true >  t4;
    _test< NT, 1, false >  t5;
    _test< NT, 2, false >  t6;
    /*_test< NT, 3, false >  t7;
    _test< NT, 4, false >  t8;
    _test< NT, 10, false > t9;*/

    void operator()() {
        t0();
        //t1();
        t2();
        t3();
        //t4();
        t5();
        t6();
        /*t7();
        t8();
        t9();*/
    }
};



int main() {
    test<unsigned int> a;
    test<int> b;
    test<float> c;
    test<double> d;
    std::cout << "-------------------------" << std::endl;
    std::cout << "type = unsigned int" << std::endl;
    std::cout << "-------------------------" << std::endl;
    a();
    std::cout << std::endl;
    std::cout << "-------------------------" << std::endl;
    std::cout << "type = signed int" << std::endl;
    std::cout << "-------------------------" << std::endl;
    b();
    std::cout << std::endl;
    std::cout << "-------------------------" << std::endl;
    std::cout << "type = float" << std::endl;
    std::cout << "-------------------------" << std::endl;
    c();
    std::cout << std::endl;
    std::cout << "-------------------------" << std::endl;
    std::cout << "type = double" << std::endl;
    std::cout << "-------------------------" << std::endl;
    d();

    if( failed != 0 ) {
        std::cout << "a total number of " << failed
                  << " tests failed!" << std::endl;
        return 1;
    }
    std::cout << "all tests passed." << std::endl;
    return 0;
}

