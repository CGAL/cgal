// enable invariant checking
#define SEGMENT_TREE_CHECK_INVARIANTS 1
#include <CGAL/Box_intersection_d.h>
// compare segment tree against brute force and simple implementations
#include <CGAL/Box_intersection_d/one_way_scan.h>
#include <CGAL/Box_intersection_d/all_pairs.h>

#include "Timer.h"

#include <iostream>
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
typedef CGAL::Default_box_d< Number_type, DIM >  Box;
typedef CGAL::Default_box_traits_d< Box > Box_adapter;
typedef CGAL::Default_box_predicate_traits_d< Box_adapter, CLOSED > Traits;
typedef std::vector< Box >     Box_container;
typedef std::pair< Box, Box >  Box_pair;
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
        //storage.push_back( std::make_pair( a, b ) );
        //std::cout << Traits::get_id( a ) << " " << Traits::get_id( b ) << std::endl;
    }
};

struct Counter_callback {
    unsigned int counter;
    Counter_callback() : counter( 0 ) {}
    void operator()( const Box& a, const Box& b ) {
        assert_intersection( a, b );
        ++counter;
        //std::cout << Traits::get_id( a ) << " " << Traits::get_id( b ) << std::endl;
    }
};

/*bool
operator==( const Box& a, const Box& b ) {
    for( unsigned int dim = 0; dim < DIM; ++dim )
        if( Traits::get_lo( a, dim ) != Traits::get_lo( b, dim ) ||
            Traits::get_hi( a, dim ) != Traits::get_hi( b, dim )   )
            return false;
    return true;
}

bool
operator==( const BoxPair& a, const BoxPair& b ) {
    return( a.first == b.first && a.second == b.second ||
            a.first == b.second && a.second == b.first );
}

template< class Storage >
unsigned int countMissingItems( Storage& a, Storage& b ) {
    unsigned int missing = 0;
    typedef typename Storage::iterator IT;
    for( IT it = a.begin(); it != a.end(); ++it ) {
        bool found = false;
        for( IT it2 = b.begin(); it2 != b.end(); ++it2 ) {
            if( *it == *it2 ) {
                found = true;
                break;
            }
        }
        if (!found) ++missing;
    }
    return missing;
}

template< class Storage >
unsigned int countDuplicates( Storage& storage ) {
    unsigned int counter = 0;
    typedef typename Storage::iterator IT;
    for( IT it = storage.begin(); it != storage.end(); ++it )
        for( IT it2 = it; it2 != storage.end(); ++it2 )
            if( it != it2 &&  *it == *it2 )
                ++counter;

    return counter;
} */

static void
test_n( unsigned int n, std::ostream& outfile )
{
    Box_container boxes1, boxes2;
    //Result_container result_allpairs, result_scanner, result_tree;
    std::cout << "generating random box sets with size " << n << " ... " << std::flush;
    fill_boxes( n, boxes1 );
    fill_boxes( n, boxes2 );
    std::cout << std::endl;
    Counter_callback callback1, callback2;
    Timer timer, timer_scan;

    std::cout << "one way scan ... " << std::flush;
    timer_scan.start();
    CGAL::one_way_scan( boxes1.begin(), boxes1.end(),
                        boxes2.begin(), boxes2.end(), callback1, Traits(), DIM - 1 );
    CGAL::one_way_scan( boxes2.begin(), boxes2.end(),
                        boxes1.begin(), boxes1.end(), callback1, Traits(), DIM - 1 );
    timer_scan.stop();
    std::cout << "got " << callback1.counter << " intersections in "
              << timer_scan.t << " seconds."
              << std::endl;

    std::cout << "segment tree ... " << std::endl;
    unsigned int problemsize = boxes1.size() + boxes2.size();
    unsigned int cutoff = 0;
    unsigned int stepsize = 500;
    float last_time = 1e30;
    if( problemsize > 3000 )
        cutoff = 500;
    while( true )
    {
        timer.reset();
        timer.start();
        CGAL::box_intersection_d( boxes1.begin(), boxes1.end(),
                                  boxes2.begin(), boxes2.end(), callback2, cutoff );
        if( problemsize < 500 ) {
            CGAL::box_intersection_d( boxes1.begin(), boxes1.end(),
                                      boxes2.begin(), boxes2.end(), callback2, cutoff );
            CGAL::box_intersection_d( boxes1.begin(), boxes1.end(),
                                      boxes2.begin(), boxes2.end(), callback2, cutoff );

            timer.stop();
            timer.t /= 3.0;
        } else
            timer.stop();
        std::cout << "cutoff = " << cutoff << " -> t = " << timer.t << std::endl;
        if( last_time < timer.t || timer.t < 1e-4) {
            if( cutoff > 2*stepsize )
                cutoff -= 2*stepsize;
            else
                cutoff = 0;
            if( problemsize < 2000 && problemsize / stepsize > 10 || problemsize / stepsize > 1000 )
                break;
            stepsize /= 2;
        }
        cutoff += stepsize;
        last_time = timer.t;

    }
    std::cout << "optimal cutoff = " << cutoff << std::endl;
    outfile << problemsize << " " << last_time << " " << timer_scan.t << std::endl;
}

void operator()() {
    std::cout << "-------------------------" << std::endl;
    std::cout << "DIM = " << DIM << std::endl;
    std::cout << "-------------------------" << std::endl;
    std::ofstream outfile( "benchmark.data" );
    outfile << "# problemsize streamingtime scanningtime" << std::endl;
    outfile.precision(9);
    outfile << std::fixed;
    for( unsigned int n = 8; n < 1000000; n = (int)(n * 1.5)) {
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

    if( failed != 0 )
        std::cout << "a total number of " << failed << " tests failed!" << std::endl;
    else
        std::cout << "all tests passed." << std::endl;
}

