// enable invariant checking
#define SEGMENT_TREE_CHECK_INVARIANTS 1
#include <CGAL/Box_intersection_d.h>

#include <CGAL/Timer.h>

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <iterator>
#include <vector>
#include <cstdio>
#include <cmath>
#include <fstream>

static unsigned int failed = 0;

template< class NT, unsigned int DIM, bool CLOSED >
struct _test {
typedef NT Number_type;
typedef CGAL::Box_intersection_d::Box_d< Number_type, DIM >  Box;
typedef CGAL::Box_intersection_d::Box_traits_d< Box > Box_adapter;
typedef CGAL::Box_intersection_d::Box_predicate_traits_d<
                                      Box_adapter, CLOSED > Traits;
typedef std::vector< Box >      Box_container;
typedef std::pair< Box, Box >   Box_pair;
typedef std::vector< Box_pair > Result_container;


static void fill_boxes( unsigned int n, Box_container& boxes ) {
    NT maxEdgeLength = (NT) pow(n, (DIM-1.0)/DIM);

    for( unsigned int i = 0; i < n; ++i ) {
        NT lo[DIM], hi[DIM];
        for( unsigned int d = 0; d < DIM; ++d ) {
            lo[d] = (NT)(drand48() * (n - maxEdgeLength));
            hi[d] = (NT)(lo[d] + 1 + (drand48() * maxEdgeLength));
        }
        boxes.push_back( Box( &lo[0], &hi[0]) );
    }
}

static void assert_intersection( const Box& a, const Box& b ) {
    for( unsigned int dim = 0; dim < DIM; ++dim ) {
        if( Traits::does_intersect( a, b, dim ) == false ) {
            std::cout << "does not intersect!" << std::endl;
            //cout << a << endl << b << endl;
            exit(-1);
        }
    }
}


template< class Storage >
struct Storage_callback {
    unsigned int counter;
    Storage& storage;
    Storage_callback( Storage& storage ) : counter( 0 ), storage( storage ) {}
    void operator()( const Box& a, const Box& b ) {
        assert_intersection( a, b );
        ++counter;
    }
};

struct Counter_callback {
    unsigned int counter;
    Counter_callback() : counter( 0 ) {}
    void operator()( const Box& a, const Box& b ) {
        assert_intersection( a, b );
        ++counter;
    }
};

static void
test_n( unsigned int n, std::ostream& outfile )
{
    Box_container boxes1, boxes2;
    //Result_container result_allpairs, result_scanner, result_tree;
    std::cout << "generating random box sets with size "
              << n << " ... " << std::flush;
    fill_boxes( n, boxes1 );
    fill_boxes( n, boxes2 );
    std::cout << std::endl;
    Counter_callback callback1, callback2;
    CGAL::Timer timer, timer_scan;
    double time, time_scan;
    unsigned int problemsize = boxes1.size() + boxes2.size();
    unsigned int cutoff = 0;
    unsigned int stepsize = 500;
    double last_time = 1e30;

    std::cout << "one way scan ... " << std::flush;
    timer_scan.start();
    unsigned int repetitions = 1;
    if( problemsize < 20 )
        repetitions = 10000;
    else if( problemsize < 50 )
        repetitions = 5000;
    else if( problemsize < 90 )
        repetitions = 2000;
    else if( problemsize < 300 )
        repetitions = 1000;
    else if( problemsize < 1000 )
        repetitions = 300;
    else if( problemsize < 10000 )
        repetitions = 30;
    else
        repetitions = 1;

    for( unsigned int i = repetitions; i; --i ) {
        CGAL::Box_intersection_d::one_way_scan( boxes1.begin(), boxes1.end(),
                                                boxes2.begin(), boxes2.end(),
                                                callback1, Traits(), DIM - 1 );
        CGAL::Box_intersection_d::one_way_scan( boxes2.begin(), boxes2.end(),
                                                boxes1.begin(), boxes1.end(),
                                                callback1, Traits(), DIM - 1 );
    }
    timer_scan.stop();
    time_scan = timer_scan.time() / repetitions;

    std::cout << "got " << callback1.counter/repetitions << " intersections in "
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
                                      CGAL::Box_intersection_d::BIPARTITE,
                                      CGAL::Box_intersection_d::CLOSED);
        }
        timer.stop();
        time = timer.time() / repetitions;
        std::cout << "cutoff = " << cutoff << " -> t = " << time << std::endl;
        if( last_time < time || time < 1e-4) {
            if( cutoff > stepsize )
                cutoff -= stepsize;
            else
                cutoff = 0;
            if( problemsize < 2000 && problemsize/stepsize > 10 ||
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



int main( int argc, char ** argv ) {
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

