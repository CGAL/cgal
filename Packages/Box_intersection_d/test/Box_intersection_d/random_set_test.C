#include <CGAL/Box_intersection_d/box_traits.h>
#include <CGAL/Box_intersection_d/one_way_scan.h>
#include <CGAL/Box_intersection_d/all_pairs.h>
// enable invariant checking
#define SEGMENT_TREE_CHECK_INVARIANTS 1
#include <CGAL/Box_intersection_d/segment_tree.h>

#include "Timer.h"

#include <iostream>
#include <cstdlib>
#include <iterator>
#include <vector>
#include <cstdio>
#include <cmath>

static unsigned int failed = 0;

template< class NT, unsigned int DIM, bool CLOSED >
struct _test {
typedef NT NumberType;
typedef CGAL::Default_Bbox_d< NumberType, DIM >  Box;
typedef CGAL::Default_Bbox_d_Adapter< Box >    BoxAdapter;
typedef CGAL::Default_Box_Traits< BoxAdapter, CLOSED > Traits;
typedef std::vector< Box >     BoxContainer;
typedef std::pair< Box, Box >  BoxPair;
typedef std::vector< BoxPair > ResultContainer;


static void fill_boxes( unsigned int n, BoxContainer& boxes ) {
    NT maxEdgeLength = (NT) pow(n, (DIM-1.0)/DIM);

    for( unsigned int i = 0; i < n; ++i ) {
        NumberType lo[DIM], hi[DIM];
        for( unsigned int d = 0; d < DIM; ++d ) {
            lo[d] = (NumberType)(drand48() * (n - maxEdgeLength));
            hi[d] = (NumberType)(lo[d] + 1 + (drand48() * maxEdgeLength));
        }
        boxes.push_back( Box( &lo[0], &hi[0]) );
    }
}

static void assertIntersection( const Box& a, const Box& b ) {
    for( unsigned int dim = 0; dim < DIM; ++dim ) {
        if( Traits::does_intersect( a, b, dim ) == false ) {
            std::cout << "does not intersect!" << std::endl;
            //cout << a << endl << b << endl;
            exit(-1);
        }
    }
}


template< class Storage >
struct StorageCallback {
    unsigned int counter;
    Storage& storage;
    StorageCallback( Storage& storage ) : counter( 0 ), storage( storage ) {}
    void operator()( const Box& a, const Box& b ) {
        assertIntersection( a, b );
        ++counter;
        //storage.push_back( std::make_pair( a, b ) );
        //std::cout << Traits::get_num( a ) << " " << Traits::get_num( b ) << std::endl;
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
test_n( unsigned int n, bool bipartite = true )
{
    BoxContainer boxes1, boxes2;
    ResultContainer result_scanner, result_tree;
    std::cout << "generating random box sets with size " << n << " ... " << std::flush;
    fill_boxes( n, boxes1 );
    if( bipartite )
        fill_boxes( n, boxes2 );
    else
        boxes2 = boxes1;
    std::cout << std::endl;
    StorageCallback< ResultContainer >
        callback1( result_scanner ),
        callback2( result_tree );

    std::cout << "one way scan ... " << std::flush;
    Timer timer;
    timer.start();
    CGAL::one_way_scan( boxes1.begin(), boxes1.end(),
                        boxes2.begin(), boxes2.end(), callback1, Traits(), DIM - 1 );
    if( bipartite )
        CGAL::one_way_scan( boxes2.begin(), boxes2.end(),
                            boxes1.begin(), boxes1.end(), callback1, Traits(), DIM - 1 );
    timer.stop();
    std::cout << "got " << callback1.counter << " intersections in "
         << timer.t << " seconds."
         << std::endl;

    std::cout << "segment tree ... " << std::flush;
    timer.reset();
    timer.start();
    Traits::cutoff = n < 200 ? 6 : n < 2000 ? 20 : n / 50;
    //Traits::cutoff = 5;
    CGAL::segment_tree( boxes1.begin(), boxes1.end(),
                        boxes2.begin(), boxes2.end(), callback2, Traits(), bipartite );
    timer.stop();
    std::cout << "got " << callback2.counter << " intersections in "
              << timer.t << " seconds." << std::endl;

    if( callback1.counter != callback2.counter ) {
        ++failed;
        /*unsigned int missing    = countMissingItems( result_scanner,
                                                     result_tree );
        unsigned int duplicates = countDuplicates( result_tree );*/
        std::cout << "!! failed !! " << std::endl;
        /*std::cout << "!! failed !! " << missing  << " missing and "
             << duplicates << " duplicate intersections in tree result. "
             << std::endl;*/
    }
    else
        std::cout << "--- passed --- " << std::endl;
}

void operator()() {
    std::cout << "-------------------------" << std::endl;
    std::cout << "DIM = " << DIM << std::endl;
    std::cout << "-------------------------" << std::endl;
    for( unsigned int n = 8; n < 200000; n = (int)(n * 1.3)) {
        std::cout << "bipartite case: " << std::endl;
        test_n( n, true );
        std::cout << "complete case: " << std::endl;
        test_n( n, false );
    }
}

}; // end struct  test

template< class NT >
struct test
{
    _test< NT, 1, true >   t0;
    _test< NT, 2, true >   t1;
    _test< NT, 3, true >   t2;
    _test< NT, 4, true >   t3;
    _test< NT, 10, true >  t4;
    _test< NT, 1, false >  t5;
    _test< NT, 2, false >  t6;
    _test< NT, 3, false >  t7;
    _test< NT, 4, false >  t8;
    _test< NT, 10, false > t9;

    void operator()() {
        t0();
        t1();
        t2();
        t3();
        t4();
        t5();
        t6();
        t7();
        t8();
        t9();
    }
};



int main( int argc, char ** argv ) {
    test<unsigned int> a;
    test<int> b;
    test<float> c;
    test<double> d;
    std::cout << "-------------------------" << std::endl;
    std::cout << "unsigned int" << std::endl;
    std::cout << "-------------------------" << std::endl;
    a();
    std::cout << std::endl;
    std::cout << "-------------------------" << std::endl;
    std::cout << "signed int" << std::endl;
    std::cout << "-------------------------" << std::endl;
    b();
    std::cout << std::endl;
    std::cout << "-------------------------" << std::endl;
    std::cout << "float" << std::endl;
    std::cout << "-------------------------" << std::endl;
    c();
    std::cout << std::endl;
    std::cout << "-------------------------" << std::endl;
    std::cout << "double" << std::endl;
    std::cout << "-------------------------" << std::endl;
    d();

    if( failed != 0 )
        std::cout << "a total number of " << failed << " tests failed!" << std::endl;
    else
        std::cout << "all tests passed." << std::endl;
}

