// enable invariant checking
#define SEGMENT_TREE_CHECK_INVARIANTS 1
#include <CGAL/Box_intersection_d.h>
// compare segment tree against brute force and simple implementations
#include <CGAL/Box_intersection_d/one_way_scan.h>
#include <CGAL/Box_intersection_d/all_pairs.h>
#include <CGAL/Timer.h>

#include <iostream>
#include <cstdlib>
#include <iterator>
#include <vector>
#include <cstdio>
#include <cmath>

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
test_n( unsigned int n,
        CGAL::Box_intersection_d::Setting
                             setting = CGAL::Box_intersection_d::BIPARTITE )
{
    Box_container boxes1, boxes2;
    //Result_container result_allpairs, result_scanner, result_tree;
    std::cout << "generating random box sets with size " << n
              << " ... " << std::flush;
    fill_boxes( n, boxes1 );
    bool bipartite = setting == CGAL::Box_intersection_d::BIPARTITE;

    if( bipartite )
        fill_boxes( n, boxes2 );
    else
        boxes2 = boxes1;
    std::cout << std::endl;
    Counter_callback callback0, callback1, callback2;
    CGAL::Timer timer;

    if( n < 5000 ) {
        std::cout << "all pairs ... " << std::flush;
        timer.start();
        CGAL::Box_intersection_d::all_pairs( boxes1.begin(), boxes1.end(),
                                             boxes2.begin(), boxes2.end(),
                                             callback0, Traits(), DIM - 1 );
        timer.stop();
        std::cout << "got " << callback0.counter << " intersections in "
                  << timer.time() << " seconds."
            << std::endl;
        timer.reset();
    }
    std::cout << "one way scan ... " << std::flush;
    timer.start();
    CGAL::Box_intersection_d::one_way_scan( boxes1.begin(), boxes1.end(),
                                            boxes2.begin(), boxes2.end(),
                                            callback1, Traits(), DIM - 1 );
    if( bipartite )
        CGAL::Box_intersection_d::one_way_scan( boxes2.begin(), boxes2.end(),
                                                boxes1.begin(), boxes1.end(),
                                                callback1, Traits(), DIM - 1);
    timer.stop();
    std::cout << "got " << callback1.counter << " intersections in "
              << timer.time() << " seconds."
              << std::endl;

    std::cout << "segment tree ... " << std::flush;
    timer.reset();
    timer.start();
    unsigned int cutoff = n < 200 ? 6 : n < 2000 ? 20 : n / 50;
    CGAL::box_intersection_d_custom( boxes1.begin(), boxes1.end(),
                                     boxes2.begin(), boxes2.end(),
                                     callback2, Traits(), cutoff, setting );
    timer.stop();
    std::cout << "got " << callback2.counter << " intersections in "
              << timer.time() << " seconds." << std::endl;

    if( callback1.counter != callback2.counter ||
        n < 20000 && callback0.counter != callback1.counter )
    {
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
    for( unsigned int n = 1024; n < 200000; n = (int)(n * 8)) {
        std::cout << "bipartite case: " << std::endl;
        test_n( n, CGAL::Box_intersection_d::BIPARTITE );
        //std::cout << "complete case: " << std::endl;
        //test_n( n, CGAL::Box_intersection_d::COMPLETE );
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
    //test<unsigned int> a;
    //test<int> b;
    test<float> c;
    //test<double> d;
    std::cout << "-------------------------" << std::endl;
    std::cout << "unsigned int" << std::endl;
    std::cout << "-------------------------" << std::endl;
    //a();
    std::cout << std::endl;
    std::cout << "-------------------------" << std::endl;
    std::cout << "signed int" << std::endl;
    std::cout << "-------------------------" << std::endl;
    //b();
    std::cout << std::endl;
    std::cout << "-------------------------" << std::endl;
    std::cout << "float" << std::endl;
    std::cout << "-------------------------" << std::endl;
    c();
    std::cout << std::endl;
    std::cout << "-------------------------" << std::endl;
    std::cout << "double" << std::endl;
    std::cout << "-------------------------" << std::endl;
    //d();

    if( failed != 0 ) {
        std::cout << "a total number of " << failed
                  << " tests failed!" << std::endl;
        return 1;
    }
    std::cout << "all tests passed." << std::endl;
    return 0;
}

