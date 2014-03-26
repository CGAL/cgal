// Copyright (c) 2004  Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
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
#include <fstream>
#include <cstdlib>

#include "util.h"

static unsigned int failed = 0;

template< class NT, unsigned int DIM, bool CLOSED >
struct _test {
    typedef Util< NT, DIM, CLOSED > Uti1;



static void
test_n( unsigned int n, std::ostream& outfile )
{
    typename Uti1::Box_container boxes1, boxes2;
    //Result_container result_allpairs, result_scanner, result_tree;
    std::cout << "generating random box sets with size "
              << n << " ... " << std::flush;
    Uti1::fill_boxes( n, boxes1 );
    Uti1::fill_boxes( n, boxes2 );
    std::cout << std::endl;
    unsigned int c1 = 0, c2 = 0;
    typename Uti1::Counter_callback callback1(c1), callback2(c2);
    CGAL::Timer timer, timer_scan;
    double time, time_scan;
    std::size_t problemsize = boxes1.size() + boxes2.size();
    unsigned int cutoff = 0;
    unsigned int stepsize = 500;
    double last_time = 1e30;

    std::cout << "one way scan ... " << std::flush;
    timer_scan.start();
    unsigned int repetitions = 1;
    if( problemsize < 20 )
        repetitions = 1000;
    else if( problemsize < 50 )
        repetitions = 500;
    else if( problemsize < 90 )
        repetitions = 100;
    else if( problemsize < 300 )
        repetitions = 50;
    else if( problemsize < 1000 )
        repetitions = 10;
    else if( problemsize < 10000 )
        repetitions = 4;
    else
        repetitions = 1;

    for( unsigned int i = repetitions; i; --i ) {
        CGAL::Box_intersection_d::one_way_scan( boxes1.begin(), boxes1.end(),
                                                boxes2.begin(), boxes2.end(),
                                                callback1, 
                                                typename Uti1::Traits(),
                                                DIM - 1 );
        CGAL::Box_intersection_d::one_way_scan( boxes2.begin(), boxes2.end(),
                                                boxes1.begin(), boxes1.end(),
                                                callback1, 
                                                typename Uti1::Traits(),
                                                DIM - 1 );
    }
    timer_scan.stop();
    time_scan = timer_scan.time() / repetitions;
    const unsigned int scan_counter = callback1.get_counter()/repetitions;

    std::cout << "got " << scan_counter << " intersections in "
              << time_scan << " seconds."
              << std::endl;

    std::cout << "segment tree ... " << std::endl;
    if( problemsize > 3000 )
        cutoff = 500;
    while( true )
    {
        timer.reset();
        timer.start();
        for( unsigned int i = repetitions; i; --i ) {
            CGAL::box_intersection_d( boxes1.begin(), boxes1.end(),
                                      boxes2.begin(), boxes2.end(),
                                      callback2, cutoff,
                                      CGAL::Box_intersection_d::CLOSED,
                                      CGAL::Box_intersection_d::BIPARTITE);
         }
        timer.stop();
        time = timer.time() / repetitions;
        const unsigned int stream_counter = callback2.get_counter()/repetitions;
        callback2.reset_counter();
        if( stream_counter != scan_counter ) {
            std::cout << "!! different number of intersections! got "
                      << stream_counter << "intersections!" << std::endl;
            std::exit(-1);
        }
        std::cout << "cutoff = " << cutoff << " -> t = " << time << std::endl;
        if( last_time < time || time < 1e-4) {
            if( cutoff > stepsize )
                cutoff -= stepsize;
            else
                cutoff = 0;
            if( ( problemsize < 2000 && problemsize/stepsize > 10 ) ||
                                        problemsize/stepsize > 50 )
                break;
            if( cutoff > stepsize )
                cutoff -= stepsize;
            else
                cutoff = 0;
            stepsize /= 2;
            last_time = 1e30;
        } else {
            last_time = time;
            //stepsize *= 2;
        }
        cutoff += stepsize;
    }

    std::cout << "optimal cutoff = " << cutoff << std::endl;
    outfile << problemsize << " " << last_time << " " << time_scan << std::endl;
}

void operator()() {
    std::cout << "-------------------------" << std::endl;
    std::cout << "DIM = " << DIM << std::endl;
    std::cout << "-------------------------" << std::endl;
    std::ofstream outfile( "benchmark.data" );
    outfile << "# problemsize streamingtime scanningtime" << std::endl;
    outfile.precision(9);
    // correct for >= g++3.0 , but g++2.95 does not conform to the standard
    //outfile << std::fixed; 
    // workaround for g++2.95: (does not work for >= 3.0)
    //outfile << setiosflags( ios::fixed );
    for( unsigned int n = 1024; n < 200000; n = (int)(n * 4)) {
        test_n( n, outfile );
    }
}

}; // end struct  test

template< class NT >
struct test
{
    //_test< NT, 1, true >   t0;
    //_test< NT, 2, true >   t1;
    _test< NT, 3, true >   t2;
    /*_test< NT, 4, true >   t3;
    _test< NT, 10, true >  t4;
    _test< NT, 1, false >  t5;
    _test< NT, 2, false >  t6;
    _test< NT, 3, false >  t7;
    _test< NT, 4, false >  t8;
    _test< NT, 10, false > t9;*/

    void operator()() {
        //t0();
        //t1();
        t2();
        /*t3();
        t4();
        t5();
        t6();
        t7();
        t8();
        t9();*/
    }
};



int main() {
    test<float> c;
    c();

    if( failed != 0 ) {
        std::cout << "a total number of " << failed << " tests failed!" 
                  << std::endl;
        return 1;
    }
    std::cout << "all tests passed." << std::endl;
    return 0;
}

